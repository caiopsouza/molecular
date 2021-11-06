mod checks;

use std::fmt::{Debug, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::mem::swap;
use std::ops::{Index, IndexMut};

type RangePos = [(f64, f64, f64)];

type VecPos = Vec<(f64, f64, f64)>;

#[derive(Default, Clone)]
pub struct Torsion([f64; 12]);

impl Torsion {
    // Assumes the rest of the matrix is already zeroed.
    fn eye() -> Self {
        let mut res = Self::default();

        res[(0, 0)] = 1f64;
        res[(1, 1)] = 1f64;
        res[(2, 2)] = 1f64;

        res
    }

    fn product_line(&mut self, one: &Self, other: &Self, line: usize) {
        self[(line, 0)] = one[(line, 0)] * other[(0, 0)] + one[(line, 1)] * other[(1, 0)] + one[(line, 2)] * other[(2, 0)];
        self[(line, 1)] = one[(line, 0)] * other[(0, 1)] + one[(line, 1)] * other[(1, 1)] + one[(line, 2)] * other[(2, 1)];
        self[(line, 2)] = one[(line, 0)] * other[(0, 2)] + one[(line, 1)] * other[(1, 2)] + one[(line, 2)] * other[(2, 2)];
        self[(line, 3)] = one[(line, 0)] * other[(0, 3)] + one[(line, 1)] * other[(1, 3)] + one[(line, 2)] * other[(2, 3)] + one[(line, 3)];
    }

    fn product(&self, other: &Self) -> Self {
        let mut res = Self::default();

        res.product_line(self, other, 0);
        res.product_line(self, other, 1);
        res.product_line(self, other, 2);

        res
    }
}

impl Debug for Torsion {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "[[{}, {}, {}, {}]\n [{}, {}, {}, {}]\n [{}, {}, {}, {}]\n [0, 0, 0, 1]]",
               self[(0, 0)],
               self[(0, 1)],
               self[(0, 2)],
               self[(0, 3)],
               self[(1, 0)],
               self[(1, 1)],
               self[(1, 2)],
               self[(1, 3)],
               self[(2, 0)],
               self[(2, 1)],
               self[(2, 2)],
               self[(2, 3)],
        )
    }
}

impl Index<(usize, usize)> for Torsion {
    type Output = f64;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.0[4 * index.0 + index.1]
    }
}

impl IndexMut<(usize, usize)> for Torsion {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.0[4 * index.0 + index.1]
    }
}

#[derive(Debug)]
pub struct Problem {
    node_count: usize,
    edge_count: usize,
    node_capacity: usize,
    edges: Vec<Vec<usize>>,
    dists: Vec<f64>,
    torsions: Vec<(Torsion, Torsion)>,
}

impl Problem {
    fn from_file(file: &str) -> Self {
        let mut data = Vec::<(usize, usize, f64)>::with_capacity(1000);

        let file = File::open(file).unwrap_or_else(|_| panic!("file: {}", file));
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
            node_count,
            edge_count,
            node_capacity,
            edges: vec![Vec::<usize>::with_capacity(20); node_capacity],
            dists: vec![0f64; node_capacity * node_capacity],
            torsions: Vec::<(Torsion, Torsion)>::with_capacity(node_capacity),
        };

        // Populate the data dict
        for (node0, node1, dist) in data {
            res.edges[node1].push(node0);
            res.set_dist(node0, node1, dist);
        }

        res.compute_torsions();

        res
    }

    fn get_dist(&self, index0: usize, index1: usize) -> f64 {
        let dist = self.dists[index0 * self.node_capacity + index1];

        if cfg!(debug_assertions) && dist.abs() < 1e-10 {
            panic!("dist between {} and {} is too close to zero", index0, index1);
        }

        dist
    }

    fn set_dist(&mut self, index0: usize, index1: usize, dist: f64) {
        self.dists[index0 * self.node_capacity + index1] = dist;
    }

    fn get_torsion(&self, index: usize) -> &(Torsion, Torsion) {
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
    fn trig_theta_functions(&self, node: usize) -> (f64, f64) {
        let d_ab = self.get_dist(node - 2, node - 1);
        let d_bi = self.get_dist(node - 2, node);
        let d_ai = self.get_dist(node - 1, node);

        // Theta
        let mut cos_theta = d_ab * d_ab + d_ai * d_ai - d_bi * d_bi;
        cos_theta /= 2f64 * d_ab * d_ai;

        // Result
        (Self::sin_from_cos(cos_theta), cos_theta)
    }

    // Returns the sin(theta), cos(theta), sin(omega) and cos(omega) of the node.
    #[allow(clippy::many_single_char_names)]
    fn trig_functions(&self, node: usize) -> (f64, f64, f64, f64) {
        let d_cb = self.get_dist(node - 3, node - 2);
        let d_ca = self.get_dist(node - 3, node - 1);
        let d_ci = self.get_dist(node - 3, node);

        let d_ab = self.get_dist(node - 2, node - 1);
        let d_bi = self.get_dist(node - 2, node);

        let d_ai = self.get_dist(node - 1, node);

        let d_ab_squared = d_ab * d_ab;
        let d_ai_squared = d_ai * d_ai;
        let d_bi_squared = d_bi * d_bi;
        let d_ca_squared = d_ca * d_ca;
        let d_cb_squared = d_cb * d_cb;
        let d_ci_squared = d_ci * d_ci;

        // Theta
        let mut cos_theta = d_ab_squared + d_ai_squared - d_bi_squared;
        cos_theta /= 2f64 * d_ab * d_ai;

        // Omega
        let a = 2f64 * d_ab_squared * (d_cb_squared + d_bi_squared - d_ci_squared);
        let b = d_cb_squared + d_ab_squared - d_ca_squared;
        let c = d_ab_squared + d_bi_squared - d_ai_squared;

        let d = 4f64 * d_cb_squared * d_ab_squared - b * b;
        let e = 4f64 * d_ab_squared * d_bi_squared - c * c;

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

fn compute_position_and_error(problem: &Problem, node: usize, positions: &RangePos, torsion: &Torsion, sign: bool) -> ((f64, f64, f64), Torsion, f64) {
    let next_torsion = problem.get_torsion(node);
    let next_torsion = if sign { &next_torsion.1 } else { &next_torsion.0 };

    let cumulative_torsion = torsion.product(next_torsion);

    let pos = (cumulative_torsion[(0, 3)], cumulative_torsion[(1, 3)], cumulative_torsion[(2, 3)]);

    let mut total_err = 0f64;
    for &neighbor in &problem.edges[node] {
        let dist = problem.get_dist(neighbor, node);
        let neighbor = positions[neighbor];

        let err = (neighbor.0 - pos.0, neighbor.1 - pos.1, neighbor.2 - pos.2);
        let err = err.0 * err.0 + err.1 * err.1 + err.2 * err.2;
        let err = err.sqrt() - dist;

        total_err += err.abs();
    }

    (pos, cumulative_torsion, total_err)
}

fn step(problem: &Problem,
        node: usize,
        positions: &mut VecPos,
        torsion: &Torsion,
        sign: bool,
        parent_error: f64,
        mut error_so_far: f64) -> (f64, (f64, f64, f64)) {
    let (pos, cumulative_torsion, total_err) = compute_position_and_error(problem, node, positions, torsion, sign);

    // End the recursion if the error is bigger than the error so far
    if parent_error + total_err >= error_so_far {
        return (parent_error + total_err, pos);
    }

    positions[node] = pos;

    // End the recursion if it's the last node
    if node >= problem.node_count {
        return (parent_error + total_err, pos);
    }

    // Start recursion
    let (err_child, pos_child) = step(problem, node + 1, positions, &cumulative_torsion, true, parent_error + total_err, error_so_far);
    if err_child < error_so_far {
        error_so_far = err_child;
        positions[node + 1] = pos_child;
    }

    let (err_child, pos_child) = step(problem, node + 1, positions, &cumulative_torsion, false, parent_error + total_err, error_so_far);
    if err_child < error_so_far {
        error_so_far = err_child;
        positions[node + 1] = pos_child;
    }

    (error_so_far, pos)
}

pub fn load_problem(problem: &'static str) -> Problem {
    let problem = "/home/caio/molecular/src/data_rem/".to_owned() + problem;
    Problem::from_file(&problem)
}

pub fn solve_first_three(problem: &Problem) -> (Torsion, VecPos) {
    let mut positions = vec![(0f64, 0f64, 0f64); problem.node_capacity];

    // Position 2 (position 1 is all zeros)
    let r2 = problem.get_dist(1, 2);
    positions[2].0 = -r2;

    // Position 3
    let r3 = problem.get_dist(2, 3);
    let (sin_t3, cos_t3) = problem.trig_theta_functions(3);
    positions[3] = (r3 * cos_t3 - r2, r3 * sin_t3, 0f64);

    // Torsion 3
    let cumulative_torsion_3 = problem.get_torsion(2).0.product(&problem.get_torsion(3).0);

    (cumulative_torsion_3, positions)
}

pub fn solve(problem: &Problem, expected_error: f64) -> (f64, VecPos) {
    let (cumulative_torsion_3, mut positions) = solve_first_three(problem);
    let (err, _) = step(problem, 4, &mut positions, &cumulative_torsion_3, true, 0f64, expected_error);
    (err, positions)
}

pub fn solve_default_error(problem: &Problem) -> (f64, VecPos) {
    solve(problem, 1e-7 * problem.edge_count as f64)
}

pub fn heuristics_greedy(problem: &Problem) -> (f64, Vec<bool>, VecPos) {
    let (mut cumulative_torsion, mut positions) = solve_first_three(problem);
    let mut solution = vec![true; problem.node_capacity];

    let mut total_err = 0f64;

    for node in 4..problem.node_capacity {
        let (pos_true, cumulative_torsion_true, err_true) = compute_position_and_error(problem, node, &positions, &cumulative_torsion, true);

        // If the error is acceptable, no need to calculate the other.
        // This will usually happen when there's no extra information about the edges and the error will be zero.
        // Tn that case, calculating the alternative will also result in zero.
        if err_true < 1e-20 {
            cumulative_torsion = cumulative_torsion_true;
            total_err += err_true;
            positions[node] = pos_true;
            solution[node] = true;
            continue;
        }

        let (pos_false, cumulative_torsion_false, err_false) = compute_position_and_error(problem, node, &positions, &cumulative_torsion, false);

        let (pos_best, err_best, sol_best) =
            if err_true <= err_false {
                cumulative_torsion = cumulative_torsion_true;
                (pos_true, err_true, true)
            } else {
                cumulative_torsion = cumulative_torsion_false;
                (pos_false, err_false, false)
            };

        total_err += err_best;
        positions[node] = pos_best;
        solution[node] = sol_best;
    }

    (total_err, solution, positions)
}

pub fn heuristic_local_search(problem: &Problem) -> (f64, Vec<bool>, VecPos) {
    let (best_err, solution, positions) = heuristics_greedy(problem);

    //for node in (5..problem.node_capacity).rev() {}

    (best_err, solution, positions)
}

pub fn solve_heuristic_error(problem: &Problem) -> (f64, VecPos) {
    let (err, _, _) = heuristic_local_search(problem);
    solve(problem, err)
}

pub fn format(problem: &Problem, positions: &RangePos) -> String {
    let mut res = Vec::<String>::with_capacity(problem.node_capacity);
    for (node, position) in positions.iter().enumerate().skip(1) {
        res.push(format!("{} {:.9} {:.9} {:.9}", node, position.0, position.1, position.2));
    }
    res.join("\n")
}

pub fn load_solve_and_format(problem: &'static str) -> String {
    let problem = load_problem(problem);
    let (err, positions) = solve_default_error(&problem);
    println!("error={:e}", err / problem.edge_count as f64);
    format(&problem, &positions)
}

#[allow(dead_code)]
fn main() {
    let actual = load_solve_and_format("1fs3.nmr");
    println!("{}", actual);
}
