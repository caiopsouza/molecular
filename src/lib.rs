use std::fmt::{Debug, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::mem::swap;
use std::ops::{Index, IndexMut};

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
        if cfg!(debug_assertions) {
            &self.0[4 * index.0 + index.1]
        } else {
            unsafe {
                self.0.get_unchecked(4 * index.0 + index.1)
            }
        }
    }
}

impl IndexMut<(usize, usize)> for Torsion {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        if cfg!(debug_assertions) {
            &mut self.0[4 * index.0 + index.1]
        } else {
            unsafe {
                self.0.get_unchecked_mut(4 * index.0 + index.1)
            }
        }
    }
}

#[derive(Debug)]
pub struct Problem {
    largest_node: usize,
    node_capacity: usize,
    edges: Vec<Vec<usize>>,
    dists: Vec<f64>,
    torsions: Vec<(Torsion, Torsion)>,
}

impl Problem {
    fn from_file(file: &str) -> Self {
        let mut data = Vec::<(usize, usize, f64)>::with_capacity(1000);

        let file = File::open(file).unwrap();
        let mut reader = BufReader::new(file);
        let mut line = String::with_capacity(200);

        let mut largest_node = 0usize;

        // Find the largest node
        while reader.read_line(&mut line).unwrap() > 0 {
            let mut line_split = line.split_whitespace();

            let mut node0 = line_split.next().unwrap().parse::<usize>().unwrap();
            let mut node1 = line_split.next().unwrap().parse::<usize>().unwrap();
            let dist = line_split.next().unwrap().parse::<f64>().unwrap();

            if node0 > node1 {
                swap(&mut node0, &mut node1);
            }

            if node1 > largest_node {
                largest_node = node1;
            }

            data.push((node0, node1, dist));

            line.clear();
        }

        let node_capacity = largest_node + 1;

        let mut res = Self {
            largest_node,
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
        if cfg!(debug_assertions) {
            self.dists[index0 * self.node_capacity + index1]
        } else {
            unsafe {
                *self.dists.get_unchecked(index0 * self.node_capacity + index1)
            }
        }
    }

    fn set_dist(&mut self, index0: usize, index1: usize, dist: f64) {
        if cfg!(debug_assertions) {
            self.dists[index0 * self.node_capacity + index1] = dist;
        } else {
            unsafe {
                *self.dists.get_unchecked_mut(index0 * self.node_capacity + index1) = dist;
            }
        }
    }

    fn get_torsion(&self, index: usize) -> &(Torsion, Torsion) {
        if cfg!(debug_assertions) {
            &self.torsions[index]
        } else {
            unsafe {
                self.torsions.get_unchecked(index)
            }
        }
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

fn step(problem: &Problem,
        node: usize,
        positions: &mut Vec<(f64, f64, f64)>,
        torsion: &Torsion,
        sign: bool,
        parent_error: f64,
        mut error_so_far: f64) -> (f64, (f64, f64, f64)) {
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
        let err = err - dist * dist;

        total_err += err.abs();
    }

    // End the recursion if the error is bigger than the error so far
    if parent_error + total_err > error_so_far {
        return (parent_error + total_err, pos);
    }

    // End the recursion if it's the last node
    if node >= problem.largest_node {
        return (parent_error + total_err, pos);
    }

    positions[node] = pos;

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
    let problem = "C:/repos/molecular/src/data/".to_owned() + problem;
    Problem::from_file(&problem)
}

pub fn solve(problem: &Problem) -> Vec<(f64, f64, f64)> {
    let mut positions = vec![(0f64, 0f64, 0f64); problem.node_capacity];

    // Position 2
    let r2 = problem.get_dist(1, 2);
    positions[2].0 = -r2;

    // Position 3
    let r3 = problem.get_dist(2, 3);
    let (sin_t3, cos_t3) = problem.trig_theta_functions(3);
    positions[3] = (r3 * cos_t3 - r2, r3 * sin_t3, 0f64);

    // Torsion 3
    let cumulative_torsion_3 = problem.get_torsion(2).0.product(&problem.get_torsion(3).0);

    step(problem, 4, &mut positions, &cumulative_torsion_3, true, 0f64, f64::MAX);

    positions
}

pub fn assert_solution(problem: &Problem, res: &[(f64, f64, f64)]) {
    let mut res_array = Vec::<String>::with_capacity(problem.node_capacity);
    for node in 1..problem.node_capacity {
        res_array.push(format!("{} {:.9} {:.9} {:.9}", node, res[node].0, res[node].1, res[node].2));
    }

    let actual = res_array.join("\n");

    let expected = "1 0.000000000 0.000000000 0.000000000
2 -1.435771918 0.000000000 0.000000000
3 -1.990953413 1.405032920 0.000000000
4 -3.230578578 1.574823282 -0.278184623
5 -3.846005018 2.937376857 -0.386973647
6 -3.209544596 3.658055379 -1.550661578
7 -3.149906990 4.991645267 -1.279631384
8 -2.638373793 5.861907983 -2.359386260
9 -3.730635718 6.056897287 -3.456112323
10 -3.242143785 6.367473414 -4.629146466
11 -4.155159971 6.372586681 -5.824640213
12 -5.204978526 7.469371064 -5.666163027
13 -6.454073161 7.168461227 -5.879566078
14 -7.476383865 8.133009456 -5.651641664
15 -7.394218307 9.294410569 -6.707608650
16 -7.827918112 10.394797522 -6.217508626
17 -7.990918235 11.645314376 -7.070482575
18 -9.287250178 11.562293185 -7.854172549
19 -9.229923524 11.638205070 -9.178151038
20 -10.424239961 11.733574803 -10.047611138
21 -11.176548170 13.040520980 -9.776255192
22 -12.399633799 13.000298073 -9.389291808
23 -13.217759565 14.183947004 -9.014747499
24 -13.611918270 15.063025610 -10.158544998
25 -13.340315242 14.763105929 -11.411128321
26 -13.723792585 15.769522544 -12.455326925
27 -14.858030538 15.086033988 -13.193404605
28 -14.954782674 15.409674095 -14.502738827
29 -16.099606569 15.005279981 -15.370057482
30 -17.414701938 15.574184364 -14.861870535
31 -17.359714096 16.722272932 -14.280173896
32 -18.570322117 17.375873866 -13.858877043
33 -18.893819181 17.167410719 -12.386343931
34 -18.293701579 16.113197831 -11.855782859
35 -18.673372607 15.880952505 -10.453427334
36 -20.024025156 15.186806989 -10.377593991
37 -20.825747202 15.472520405 -9.356555498
38 -22.212952910 14.893303415 -9.191889164
39 -22.135693424 13.364202615 -9.168571486
40 -23.166278423 12.669706631 -9.558828951
41 -23.071269600 11.160460243 -9.604596462
42 -22.663688150 10.559800779 -8.265372221
43 -22.954105447 11.158540978 -7.180032197
44 -22.582938556 10.501849596 -5.879322941
45 -21.117620162 10.552986944 -5.684002749
46 -20.470340473 11.620863492 -6.137633306
47 -19.011869268 11.684832260 -6.174602740
48 -18.409539616 10.790198615 -7.246952325
49 -18.988996554 10.497328773 -8.345259163
50 -18.321453201 9.475462123 -9.199553856
51 -18.539033024 8.069879191 -8.643833025
52 -19.603204834 7.797896598 -7.871690662
53 -19.776946218 6.486032674 -7.182557166
54 -18.741855629 6.422493338 -6.038901962
55 -18.517916161 7.481953957 -5.254019608
56 -17.385729370 7.361311343 -4.250312764
57 -16.066993450 7.087128287 -4.975230157
58 -15.698195319 7.766662470 -6.048402665
59 -14.443252262 7.564841692 -6.728340380
60 -14.404655603 6.117819477 -7.224564127
61 -15.498014499 5.570291646 -7.815055954
62 -15.589471226 4.243812348 -8.331516348
63 -15.202403480 3.266230164 -7.218341772
64 -15.816457139 3.462584831 -6.013995913
65 -15.528644016 2.520748195 -4.897694719
66 -14.066663896 2.636174891 -4.492370395
67 -13.475699559 3.862372826 -4.454045751
68 -12.081273345 3.960723283 -4.012144494
69 -11.213157048 3.378465113 -5.157569533
70 -11.547011607 3.564517784 -6.399857961
71 -10.625869479 2.952077375 -7.481015639
72 -10.721348432 1.409918386 -7.338501001
73 -11.879240559 0.844015431 -7.029744354
74 -12.037979560 -0.563795680 -6.847951047
75 -11.088943029 -1.080717717 -5.768854540
76 -11.095577788 -0.464773071 -4.633061789
77 -10.169962805 -0.795395709 -3.554947429
78 -8.741758245 -0.830420694 -4.056801975
79 -8.299492316 0.316054834 -4.639046522
80 -6.968470323 0.523986461 -5.146157843
81 -6.605627177 -0.552741477 -6.159364534
82 -7.393183324 -0.816543552 -7.106662154
83 -7.116066886 -1.807569500 -8.089240420
84 -7.033601141 -3.233146449 -7.429059565
85 -7.931137849 -3.554340184 -6.529500538
86 -7.730685398 -4.891412420 -5.823558861
87 -6.340842781 -4.982689355 -5.194429288
88 -5.832150564 -3.893988525 -4.517388866
89 -4.600619303 -3.919994061 -3.799759807
90 -3.464740422 -4.009656774 -4.744782963
91 -3.514393155 -3.257913733 -5.862517280
92 -2.362137717 -3.294467398 -6.823200493
93 -2.234928793 -4.634984002 -7.460447010
94 -3.360556047 -5.263165313 -7.790203781
95 -3.387737941 -6.594193830 -8.378886404
96 -3.310199859 -7.723736255 -7.293648755
97 -3.114383938 -7.444507655 -6.066938310
98 -2.992285158 -8.429369288 -4.951904476
99 -4.183633853 -9.384390861 -4.951615591
100 -5.318901913 -8.829502666 -5.221359030
101 -6.564802447 -9.696089423 -5.269900259
102 -7.684683663 -9.021791576 -4.484237813
103 -7.537648470 -9.093887240 -3.115979845
104 -8.402082427 -8.444374870 -2.198452399
105 -9.717934880 -9.210047272 -2.048789132
106 -10.092857243 -10.116087223 -2.806759336
107 -11.596166350 -10.398505800 -2.997755345
108 -12.192282620 -9.222488317 -3.767175434";

    //println!("{}", actual);
    assert!(actual == expected);
    //println!("Tests passed.");
}

fn main() {
    let problem = load_problem("1ppt.nmr");
    let positions = solve(&problem);
    assert_solution(&problem, &positions);
}
