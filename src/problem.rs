use std::fs::File;
use std::io::{BufRead, BufReader};
use std::mem::swap;
use crate::torsion::Torsion;

#[derive(Debug)]
pub struct Problem {
    pub(crate) data: Vec<(usize, usize, f64)>,
    pub node_count: usize,
    pub(crate) edge_count: usize,
    pub(crate) node_capacity: usize,
    pub edges: Vec<Vec<usize>>,
    dists: Vec<f64>,
    torsions: Vec<(Torsion, Torsion)>,
}

impl Problem {
    pub(crate) fn from_file(file: &str) -> Self {
        let mut data = Vec::<(usize, usize, f64)>::with_capacity(1000);

        let file = File::open(file).unwrap_or_else(|_| panic!("file: {} not found", file));
        let mut reader = BufReader::new(file);
        let mut line = String::with_capacity(200);

        let mut node_count = 0usize;
        let mut edge_count = 0usize;

        // Find the largest node and the edge count
        while reader.read_line(&mut line).unwrap() > 0 {
            edge_count += 1;

            let mut line_split = line.split_whitespace();

            let mut node0 = line_split.next().unwrap().parse::<usize>().unwrap();
            let mut node1 = line_split.next().unwrap().parse::<usize>().unwrap();
            let dist = line_split.next().unwrap().parse::<f64>().unwrap();

            if node0 > node1 {
                swap(&mut node0, &mut node1);
            }

            if node1 > node_count {
                node_count = node1;
            }

            data.push((node0, node1, dist));

            line.clear();
        }

        let node_capacity = node_count + 1;

        let mut res = Self {
            data: Vec::new(),
            node_count,
            edge_count,
            node_capacity,
            edges: vec![Vec::<usize>::with_capacity(20); node_capacity],
            dists: vec![0f64; node_capacity * node_capacity],
            torsions: Vec::<(Torsion, Torsion)>::with_capacity(node_capacity),
        };

        // Populate the instances_test dict
        for &(node0, node1, dist) in &data {
            res.set_dist(node0, node1, dist);

            if (node0 as isize) < (node1 as isize - 3) {
                res.edges[node1].push(node0);
            }
        }

        res.data = data;

        res.compute_torsions();

        res
    }

    pub fn get_dist(&self, index0: usize, index1: usize) -> f64 {
        let dist = self.dists[index0 * self.node_capacity + index1];

        if cfg!(debug_assertions) && dist.abs() < 1e-10 {
            panic!("dist between {} and {} is too close to zero", index0, index1);
        }

        dist
    }

    fn set_dist(&mut self, index0: usize, index1: usize, dist: f64) {
        self.dists[index0 * self.node_capacity + index1] = dist;
    }

    pub fn get_torsion(&self, index: usize) -> &(Torsion, Torsion) {
        &self.torsions[index]
    }

    fn compute_torsions(&mut self) {
        self.torsions.push((Torsion::default(), Torsion::default()));
        self.torsions.push((Torsion::eye(), Torsion::eye()));

        let mut t2 = Torsion::default();
        t2[(0, 0)] = -1f64;
        t2[(0, 3)] = -self.get_dist(1, 2);
        t2[(1, 1)] = 1f64;
        t2[(2, 2)] = -1f64;
        self.torsions.push((t2.clone(), t2));

        let mut t3 = Torsion::default();

        let (sin_theta_3, cos_theta_3) = self.trig_theta_functions(3);

        t3[(0, 0)] = -cos_theta_3;
        t3[(0, 1)] = -sin_theta_3;
        t3[(0, 3)] = -self.get_dist(2, 3) * cos_theta_3;
        t3[(1, 0)] = sin_theta_3;
        t3[(1, 1)] = -cos_theta_3;
        t3[(1, 3)] = self.get_dist(2, 3) * sin_theta_3;
        t3[(2, 2)] = 1f64;
        self.torsions.push((t3.clone(), t3));

        for node in 4..self.node_capacity {
            let torsion_pair = self.compute_torsion(node);
            self.torsions.push(torsion_pair);
        }
    }

    fn sin_from_cos(v: f64) -> f64 {
        (1f64 - v * v).sqrt()
    }

    // Returns the sin(theta), cos(theta) of the node.
    pub(crate) fn trig_theta_functions(&self, node: usize) -> (f64, f64) {
        let d_ba = self.get_dist(node - 2, node - 1);
        let d_bi = self.get_dist(node - 2, node);
        let d_ai = self.get_dist(node - 1, node);

        // Theta
        let mut cos_theta = d_ba * d_ba + d_ai * d_ai - d_bi * d_bi;
        cos_theta /= 2f64 * d_ba * d_ai;

        // Result
        (Self::sin_from_cos(cos_theta), cos_theta)
    }

    // Returns the sin(theta), cos(theta), sin(omega) and cos(omega) of the node.
    #[allow(clippy::many_single_char_names)]
    fn trig_functions(&self, node: usize) -> (f64, f64, f64, f64) {
        let d_cb = self.get_dist(node - 3, node - 2);
        let d_ca = self.get_dist(node - 3, node - 1);
        let d_ci = self.get_dist(node - 3, node);

        let d_ba = self.get_dist(node - 2, node - 1);
        let d_bi = self.get_dist(node - 2, node);

        let d_ai = self.get_dist(node - 1, node);

        let d_ba_squared = d_ba * d_ba;
        let d_ai_squared = d_ai * d_ai;
        let d_bi_squared = d_bi * d_bi;
        let d_ca_squared = d_ca * d_ca;
        let d_cb_squared = d_cb * d_cb;
        let d_ci_squared = d_ci * d_ci;

        // Theta
        let mut cos_theta = d_ba_squared + d_ai_squared - d_bi_squared;
        cos_theta /= 2f64 * d_ba * d_ai;

        // Omega
        let a = 2f64 * d_ba_squared * (d_cb_squared + d_bi_squared - d_ci_squared);
        let b = d_cb_squared + d_ba_squared - d_ca_squared;
        let c = d_ba_squared + d_bi_squared - d_ai_squared;

        let d = 4f64 * d_cb_squared * d_ba_squared - b * b;
        let e = 4f64 * d_ba_squared * d_bi_squared - c * c;

        let cos_omega = (a - b * c) / (d * e).sqrt();

        // Result
        (Self::sin_from_cos(cos_theta), cos_theta, Self::sin_from_cos(cos_omega), cos_omega)
    }

    fn compute_torsion(&self, node: usize) -> (Torsion, Torsion) {
        let ri = self.get_dist(node - 1, node);
        let (sin_ti, cos_ti, sin_oi, cos_oi) = self.trig_functions(node);

        let mut res = Torsion::default();

        res[(0, 0)] = -cos_ti;
        res[(0, 1)] = -sin_ti;
        res[(0, 3)] = -ri * cos_ti;

        res[(1, 0)] = sin_ti * cos_oi;
        res[(1, 1)] = -cos_ti * cos_oi;
        res[(1, 2)] = -sin_oi;
        res[(1, 3)] = ri * sin_ti * cos_oi;

        res[(2, 0)] = sin_ti * sin_oi;
        res[(2, 1)] = -cos_ti * sin_oi;
        res[(2, 2)] = cos_oi;
        res[(2, 3)] = ri * sin_ti * sin_oi;

        let mut res_neg = res.clone();
        res_neg[(1, 2)] *= -1f64;
        res_neg[(2, 0)] *= -1f64;
        res_neg[(2, 1)] *= -1f64;
        res_neg[(2, 3)] *= -1f64;

        (res_neg, res)
    }
}
