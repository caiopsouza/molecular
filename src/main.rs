mod tests;
mod torsion;
mod problem;
mod solution;

use std::cmp::Ordering;
use rand::prelude::StdRng;
use rand::{Rng, SeedableRng};
use crate::problem::Problem;
use crate::solution::{Pos, VecPos};
use crate::torsion::Torsion;

fn compute_position_and_error(problem: &Problem, node: usize, positions: &VecPos, torsion: &Torsion, sign: bool) -> (Pos, Torsion, f64) {
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

        total_err += err.abs() / dist;
    }

    (pos, cumulative_torsion, total_err)
}

fn step(problem: &Problem,
        node: usize,
        positions: &mut VecPos,
        torsion: &Torsion,
        sign: bool,
        parent_error: f64,
        mut error_so_far: f64) -> f64 {
    let (pos, cumulative_torsion, error) = compute_position_and_error(problem, node, positions, torsion, sign);
    let cumulative_error = parent_error + error;

    // End the recursion if the error is bigger than the error so far
    if cumulative_error >= error_so_far {
        return cumulative_error;
    }

    let mut pos_to_revert = positions[node];
    positions[node] = pos;

    // End the recursion if it's the last node
    if node >= problem.node_count {
        return cumulative_error;
    }

    // Start recursion
    let err_child = step(problem, node + 1, positions, &cumulative_torsion, true, cumulative_error, error_so_far);
    if err_child < error_so_far {
        error_so_far = err_child;
        pos_to_revert = pos;
    }

    let err_child = step(problem, node + 1, positions, &cumulative_torsion, false, cumulative_error, error_so_far);
    if err_child < error_so_far {
        error_so_far = err_child;
        pos_to_revert = pos;
    }

    positions[node] = pos_to_revert;

    error_so_far
}

pub fn load_problem_test(problem: &str) -> Problem {
    let problem = "/home/caio/molecular/src/instances/".to_owned() + problem + ".nmr";
    Problem::from_file(&problem)
}

pub fn load_problem(problem: &str) -> Problem {
    let problem = "/home/caio/molecular/src/instances/".to_owned() + problem + ".nmr";
    Problem::from_file(&problem)
}

pub fn solve_first_three(problem: &Problem) -> (Torsion, VecPos) {
    let mut positions = VecPos::with_size(problem.node_capacity);

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
    let total_error = step(problem, 4, &mut positions, &cumulative_torsion_3, true, 0f64, expected_error);
    (total_error, positions)
}

pub fn solve_default_error(problem: &Problem) -> (f64, VecPos) {
    solve(problem, 1e-7 * problem.edge_count as f64)
}

fn calculate_random_guess(problem: &Problem, node: usize, positions: &VecPos, cumulative_torsion: &Torsion,  random: &mut StdRng) -> (Pos, f64, Torsion, bool) {
    let guess = random.gen_bool(0.5);
    let (pos, cumulative_torsion, err) = compute_position_and_error(problem, node, &positions, &cumulative_torsion, guess);
    (pos, err, cumulative_torsion, guess)
}

fn calculate_best_guess(problem: &Problem, node: usize, positions: &VecPos, cumulative_torsion: &Torsion) -> (Pos, f64, Torsion, bool) {
    let (pos_true, cumulative_torsion_true, err_true) = compute_position_and_error(problem, node, &positions, &cumulative_torsion, true);
    let (pos_false, cumulative_torsion_false, err_false) = compute_position_and_error(problem, node, &positions, &cumulative_torsion, false);

    if err_true <= err_false {
        (pos_true, err_true, cumulative_torsion_true, true)
    } else {
        (pos_false, err_false, cumulative_torsion_false, false)
    }
}

fn calculate_worst_guess(problem: &Problem, node: usize, positions: &VecPos, cumulative_torsion: &Torsion) -> (Pos, f64, Torsion, bool) {
    let (pos_true, cumulative_torsion_true, err_true) = compute_position_and_error(problem, node, &positions, &cumulative_torsion, true);
    let (pos_false, cumulative_torsion_false, err_false) = compute_position_and_error(problem, node, &positions, &cumulative_torsion, false);

    if err_true > err_false {
        (pos_true, err_true, cumulative_torsion_true, true)
    } else {
        (pos_false, err_false, cumulative_torsion_false, false)
    }
}

pub fn random_solution(problem: &Problem, random: &mut StdRng) -> (f64, Vec<bool>, VecPos, Vec<Torsion>, Vec<f64>) {
    let (mut cumulative_torsion, mut positions) = solve_first_three(problem);
    let mut solution = vec![true; problem.node_capacity];
    let mut errors = vec![0.0; problem.node_capacity];

    let mut torsions = Vec::<Torsion>::with_capacity(problem.node_capacity);
    for node in 0..3 {
        torsions.push(problem.get_torsion(node).0.clone());
    }
    torsions.push(cumulative_torsion.clone());

    let mut total_err = 0f64;

    for node in 4..problem.node_capacity {
        let (pos_best, err_best, cumulative_torsion_guess, sol_best) =
            calculate_random_guess(problem, node, &positions, &cumulative_torsion, random);

        cumulative_torsion = cumulative_torsion_guess;
        torsions.push(cumulative_torsion.clone());
        errors[node] = err_best;
        total_err += err_best;
        positions[node] = pos_best;
        solution[node] = sol_best;
    }

    (total_err, solution, positions, torsions, errors)
}

pub fn heuristics_greedy_look_one(problem: &Problem) -> (f64, Vec<bool>, VecPos, Vec<Torsion>, Vec<f64>) {
    let (mut cumulative_torsion, mut positions) = solve_first_three(problem);
    let mut solution = vec![true; problem.node_capacity];
    let mut errors = vec![0.0; problem.node_capacity];

    let mut torsions = Vec::<Torsion>::with_capacity(problem.node_capacity);
    for node in 0..3 {
        torsions.push(problem.get_torsion(node).0.clone());
    }
    torsions.push(cumulative_torsion.clone());

    let mut total_err = 0f64;

    for node in 4..problem.node_capacity {
        let (pos_best, err_best, cumulative_torsion_guess, sol_best) =
            calculate_best_guess(problem, node, &positions, &cumulative_torsion);

        cumulative_torsion = cumulative_torsion_guess;
        torsions.push(cumulative_torsion.clone());
        errors[node] = err_best;
        total_err += err_best;
        positions[node] = pos_best;
        solution[node] = sol_best;
    }

    (total_err, solution, positions, torsions, errors)
}

pub fn heuristics_greedy_look_worst_one(problem: &Problem) -> (f64, Vec<bool>, VecPos, Vec<Torsion>, Vec<f64>) {
    let (mut cumulative_torsion, mut positions) = solve_first_three(problem);
    let mut solution = vec![true; problem.node_capacity];
    let mut errors = vec![0.0; problem.node_capacity];

    let mut torsions = Vec::<Torsion>::with_capacity(problem.node_capacity);
    for node in 0..3 {
        torsions.push(problem.get_torsion(node).0.clone());
    }
    torsions.push(cumulative_torsion.clone());

    let mut total_err = 0f64;

    for node in 4..problem.node_capacity {
        let (pos_best, err_best, cumulative_torsion_guess, sol_best) =
            calculate_worst_guess(problem, node, &positions, &cumulative_torsion);

        cumulative_torsion = cumulative_torsion_guess;
        torsions.push(cumulative_torsion.clone());
        errors[node] = err_best;
        total_err += err_best;
        positions[node] = pos_best;
        solution[node] = sol_best;
    }

    (total_err, solution, positions, torsions, errors)
}

pub fn heuristics_greedy_look_one_backup(problem: &Problem) -> (f64, Vec<bool>, VecPos, Vec<Torsion>, Vec<f64>) {
    let (mut cumulative_torsion, mut positions) = solve_first_three(problem);
    let mut solution = vec![true; problem.node_capacity];
    let mut errors = vec![0.0; problem.node_capacity];

    let mut torsions = Vec::<Torsion>::with_capacity(problem.node_capacity);
    for node in 0..3 {
        torsions.push(problem.get_torsion(node).0.clone());
    }
    torsions.push(cumulative_torsion.clone());

    let mut total_err = 0f64;

    for node in 4..problem.node_capacity {
        let (pos_best, err_best, cumulative_torsion_guess, sol_best) =
            calculate_best_guess(problem, node, &positions, &cumulative_torsion);

        cumulative_torsion = cumulative_torsion_guess;
        torsions.push(cumulative_torsion.clone());
        errors[node] = err_best;
        total_err += err_best;
        positions[node] = pos_best;
        solution[node] = sol_best;
    }

    (total_err, solution, positions, torsions, errors)
}

pub fn heuristics_greedy_look_two(problem: &Problem) -> (f64, Vec<bool>, VecPos, Vec<Torsion>, Vec<f64>) {
    let (mut cumulative_torsion, mut positions) = solve_first_three(problem);
    let mut solution = vec![true; problem.node_capacity];
    let mut errors = vec![0.0; problem.node_capacity];

    let mut torsions = Vec::<Torsion>::with_capacity(problem.node_capacity);
    for node in 0..3 {
        torsions.push(problem.get_torsion(node).0.clone());
    }
    torsions.push(cumulative_torsion.clone());

    let mut total_err = 0f64;

    for node in (4..problem.node_capacity - 1).step_by(2) {
        let (pos_true, cumulative_torsion_true, err_true) = compute_position_and_error(problem, node, &positions, &cumulative_torsion, true);
        positions[node] = pos_true;
        let (pos_true_next, err_true_next, cumulative_torsion_true_next, best_guess_true) =
            calculate_best_guess(problem, node + 1, &positions, &cumulative_torsion_true);

        let (pos_false, cumulative_torsion_false, err_false) = compute_position_and_error(problem, node, &positions, &cumulative_torsion, false);
        positions[node] = pos_false;
        let (pos_false_next, err_false_next, cumulative_torsion_false_next, best_guess_false) =
            calculate_best_guess(problem, node + 1, &positions, &cumulative_torsion_false);

        let (best_index, _) =
            [err_true + err_true_next, err_true + err_false_next, err_false + err_true_next, err_false + err_false_next]
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                .unwrap();

        if best_index == 0 || best_index == 1 {
            torsions.push(cumulative_torsion_true);
            errors[node] = err_true;
            total_err += err_true;
            positions[node] = pos_true;
            solution[node] = true;
        } else {
            torsions.push(cumulative_torsion_false);
            errors[node] = err_false;
            total_err += err_false;
            positions[node] = pos_false;
            solution[node] = false;
        }

        if best_index == 0 || best_index == 2 {
            cumulative_torsion = cumulative_torsion_true_next;
            torsions.push(cumulative_torsion.clone());
            errors[node + 1] = err_true_next;
            total_err += err_true_next;
            positions[node + 1] = pos_true_next;
            solution[node + 1] = best_guess_true;
        } else {
            cumulative_torsion = cumulative_torsion_false_next;
            torsions.push(cumulative_torsion.clone());
            errors[node + 1] = err_false_next;
            total_err += err_false_next;
            positions[node + 1] = pos_false_next;
            solution[node + 1] = best_guess_false;
        }
    }

    if problem.node_capacity % 2 == 0 {
        let (pos_best, err_best, cumulative_torsion_guess, sol_best) =
            calculate_best_guess(problem, problem.node_capacity - 1, &positions, &cumulative_torsion);

        cumulative_torsion = cumulative_torsion_guess;
        torsions.push(cumulative_torsion.clone());
        errors[problem.node_capacity - 1] = err_best;
        total_err += err_best;
        positions[problem.node_capacity - 1] = pos_best;
        solution[problem.node_capacity - 1] = sol_best;
    }

    (total_err, solution, positions, torsions, errors)
}

fn improve_backward(problem: &Problem, mut best_err: f64, solution: &mut Vec<bool>, mut positions: VecPos, torsions: &Vec<Torsion>) -> (f64, VecPos) {
    loop {
        let mut has_improved = false;

        for node in (5..problem.node_capacity).rev() {
            solution[node] = !solution[node];

            let mut positions_new = positions.clone();
            let mut torsions_new = torsions.clone();

            for node2 in node..problem.node_capacity {
                let (pos, torsion, _) = compute_position_and_error(problem, node2, &positions_new, &torsions_new[node2 - 1], solution[node2]);
                positions_new[node2] = pos;
                torsions_new[node2] = torsion;
            }

            let new_error = compute_error(problem, &positions_new);

            if new_error >= best_err {
                solution[node] = !solution[node];
                continue;
            }

            positions = positions_new;
            best_err = new_error;
            has_improved = true;
            break;
        }

        if !has_improved {
            break;
        }
    }

    (best_err, positions)
}

pub fn heuristic_local_search_backward(problem: &Problem) -> (f64, Vec<bool>, VecPos, Vec<Torsion>) {
    let (best_err, mut solution, positions, torsions, _) = heuristics_greedy_look_one(problem);
    let (best_err, positions) = improve_backward(problem, best_err, &mut solution, positions, &torsions);
    (best_err, solution, positions, torsions)
}

pub fn heuristic_local_search_forward(problem: &Problem) -> (f64, Vec<bool>, VecPos) {
    let (mut best_err, mut solution, mut positions, torsions, _) = heuristics_greedy_look_one(problem);

    loop {
        let mut has_improved = false;

        for node in 5..problem.node_capacity {
            solution[node] = !solution[node];

            let mut positions_new = positions.clone();
            let mut torsions_new = torsions.clone();

            for node2 in node..problem.node_capacity {
                let (pos, torsion, _) = compute_position_and_error(problem, node2, &positions_new, &torsions_new[node2 - 1], solution[node2]);
                positions_new[node2] = pos;
                torsions_new[node2] = torsion;
            }

            let new_error = compute_error(problem, &positions_new);

            if new_error >= best_err {
                solution[node] = !solution[node];
                continue;
            }

            positions = positions_new;
            best_err = new_error;
            has_improved = true;
            break;
        }

        if !has_improved {
            break;
        }
    }

    (best_err, solution, positions)
}

pub fn heuristic_local_search_best(problem: &Problem) -> (f64, Vec<bool>, VecPos) {
    let (mut best_err, mut solution, mut positions, torsions, _) = heuristics_greedy_look_one(problem);

    loop {
        let mut found_so_far = None;
        let mut error_so_far = best_err;

        for node in 5..problem.node_capacity {
            solution[node] = !solution[node];

            let mut positions_new = positions.clone();
            let mut torsions_new = torsions.clone();

            for node2 in node..problem.node_capacity {
                let (pos, torsion, _) = compute_position_and_error(problem, node2, &positions_new, &torsions_new[node2 - 1], solution[node2]);
                positions_new[node2] = pos;
                torsions_new[node2] = torsion;
            }

            let new_error = compute_error(problem, &positions_new);

            solution[node] = !solution[node];

            if new_error < error_so_far {
                found_so_far = Some(node);
                error_so_far = new_error;
            }
        }

        match found_so_far {
            None => { return (best_err, solution, positions); }
            Some(node) => {
                best_err = error_so_far;
                solution[node] = !solution[node];

                for node2 in node..problem.node_capacity {
                    let (pos, _, _) = compute_position_and_error(problem, node2, &positions, &torsions[node2 - 1], solution[node2]);
                    positions[node2] = pos;
                }
            }
        }
    }
}

pub fn heuristic_iterated_local_search(problem: &Problem, chance_perturbation: u32, iterations: usize, rng_seed: u64) -> (f64, Vec<bool>, VecPos, Vec<f64>) {
    let (_, mut solution, mut positions, mut torsions) = heuristic_local_search_backward(&problem);
    let mut random = StdRng::seed_from_u64(rng_seed);

    let mut errors = Vec::<f64>::with_capacity(iterations);

    let mut best_err = compute_error(problem, &positions);
    errors.push(best_err);

    for _ in 1..iterations {
        let mut solution_candidate = solution.clone();
        let mut positions_candidate = positions.clone();

        let mut cumulative_torsion = torsions[4].clone();

        for node2 in 5..problem.node_capacity {
            if random.gen_ratio(chance_perturbation, 100) {
                solution_candidate[node2] = !solution_candidate[node2];
            }

            let (pos, cumulative_torsion_new, _) = compute_position_and_error(problem, node2, &positions_candidate, &cumulative_torsion, solution_candidate[node2]);
            positions_candidate[node2] = pos;
            torsions[node2] = cumulative_torsion_new.clone();
            cumulative_torsion = cumulative_torsion_new;
        }

        let (_, positions_candidate) = improve_backward(problem, best_err, &mut solution_candidate, positions_candidate, &torsions);

        let err_candidate = compute_error(problem, &positions_candidate);

        if err_candidate < best_err {
            best_err = err_candidate;
            solution = solution_candidate;
            positions = positions_candidate;
        }

        errors.push(best_err);
    }

    (best_err, solution, positions, errors)
}


pub fn heuristic_iterated_local_search2(problem: &Problem, chance_perturbation: u32, iterations: usize, rng_seed: u64) -> (f64, Vec<bool>, VecPos, Vec<f64>) {
    let mut random = StdRng::seed_from_u64(rng_seed);

    let (_, mut solution, mut positions, mut torsions, _) = random_solution(problem, &mut random);

    let mut errors = Vec::<f64>::with_capacity(iterations);

    let mut best_err = compute_error(problem, &positions);
    errors.push(best_err);

    for _ in 1..iterations {
        let mut solution_candidate = solution.clone();
        let mut positions_candidate = positions.clone();

        let mut cumulative_torsion = torsions[4].clone();

        for node2 in 5..problem.node_capacity {
            if random.gen_ratio(chance_perturbation, 100) {
                solution_candidate[node2] = !solution_candidate[node2];
            }

            let (pos, cumulative_torsion_new, _) = compute_position_and_error(problem, node2, &positions_candidate, &cumulative_torsion, solution_candidate[node2]);
            positions_candidate[node2] = pos;
            torsions[node2] = cumulative_torsion_new.clone();
            cumulative_torsion = cumulative_torsion_new;
        }

        let (_, positions_candidate) = improve_backward(problem, best_err, &mut solution_candidate, positions_candidate, &torsions);

        let err_candidate = compute_error(problem, &positions_candidate);

        if err_candidate < best_err {
            best_err = err_candidate;
            solution = solution_candidate;
            positions = positions_candidate;
        }

        errors.push(best_err);
    }

    (best_err, solution, positions, errors)
}

pub fn format(problem: &Problem, positions: &VecPos) -> String {
    let mut res = Vec::<String>::with_capacity(problem.node_capacity);
    for (node, position) in positions.0.iter().enumerate().skip(1) {
        res.push(format!("{} {:.9} {:.9} {:.9}", node, position.0, position.1, position.2));
    }
    res.join("\n")
}

pub fn solve_greedy_one(problem: &str) -> (f64, String) {
    let problem = load_problem(problem);
    let (error, _, positions, _, _) = heuristics_greedy_look_one(&problem);
    (error / problem.edge_count as f64, format(&problem, &positions))
}

pub fn solve_greedy_two(problem: &str) -> (f64, String) {
    let problem = load_problem(problem);
    let (error, _, positions, _, _) = heuristics_greedy_look_two(&problem);
    (error / problem.edge_count as f64, format(&problem, &positions))
}

pub fn solve_greedy_worst_one(problem: &str) -> (f64, String) {
    let problem = load_problem(problem);
    let (error, _, positions, _, _) = heuristics_greedy_look_worst_one(&problem);
    (error / problem.edge_count as f64, format(&problem, &positions))
}

pub fn solve_exact(problem: &str) -> (f64, String) {
    let problem = load_problem(problem);
    let (error, positions) = solve_default_error(&problem);

    (error / problem.edge_count as f64, format(&problem, &positions))
}

pub fn solve_exact_with_greedy_error(problem: &str) -> (f64, String) {
    let problem = load_problem(problem);
    let (error_candidate, _, _, _, _) = heuristics_greedy_look_one(&problem);
    let (error, positions) = solve(&problem, error_candidate);
    (error / problem.edge_count as f64, format(&problem, &positions))
}

pub fn solve_local_search_forward(problem: &str) -> (f64, String) {
    let problem = load_problem(problem);
    let (error, _, positions) = heuristic_local_search_forward(&problem);
    (error / problem.edge_count as f64, format(&problem, &positions))
}

pub fn solve_local_search_backward(problem: &str) -> (f64, String) {
    let problem = load_problem(problem);
    let (error, _, positions, _) = heuristic_local_search_backward(&problem);
    (error / problem.edge_count as f64, format(&problem, &positions))
}

#[allow(dead_code)]
fn compute_error(problem: &Problem, positions: &VecPos) -> f64 {
    let mut total_error = 0f64;

    for &(node0, node1, dist) in &problem.data {
        let pos0 = positions[node0];
        let pos1 = positions[node1];

        let diff = (pos0.0 - pos1.0, pos0.1 - pos1.1, pos0.2 - pos1.2);
        let diff = diff.0 * diff.0 + diff.1 * diff.1 + diff.2 * diff.2;
        let diff = diff.sqrt();

        total_error += (diff - dist).abs() / dist;
    }

    total_error
}

#[allow(dead_code)]
fn main() {
    let instances = ["1ppt", "2erl", "1crn", "1jk2", "1pht", "1a70", "1fs3", "1hoe", "1poa", "1mbn", "1ptq", "1m40", "1n4w", "1mqq", "1bpm", "3b34", "2e7z", "1rwh", "1rgs"];

    println!("inst || exact        || greedy_one   | greedy_two   | greedy_worst || local_front  | local_back   | local_best   || ils 5%       | ils 10%      | ils 20%      | ils 30%      | ils 40%      | ils 50%      ");

    let mut rng_seed = 666u64;
    let iterations = 5;

    for file in instances.iter().copied() {
        let problem = load_problem(file);

        let (_, positions_exact) = solve_default_error(&problem);

        let (_, _, positions_greedy_one, _, _) = heuristics_greedy_look_one(&problem);
        let (_, _, positions_greedy_two, _, _) = heuristics_greedy_look_two(&problem);
        let (_, _, positions_greedy_worst, _, _) = heuristics_greedy_look_worst_one(&problem);

        let (_, _, positions_local_forward) = heuristic_local_search_forward(&problem);
        let (_, _, positions_local_backward, _) = heuristic_local_search_backward(&problem);
        let (_, _, positions_local_best) = heuristic_local_search_best(&problem);

        let (_, _, positions_iterated_5, _) = heuristic_iterated_local_search(&problem, 5, iterations, rng_seed);
        rng_seed += 1;

        let (_, _, positions_iterated_10, _) = heuristic_iterated_local_search(&problem, 10, iterations, rng_seed);
        rng_seed += 1;

        let (_, _, positions_iterated_20, _) = heuristic_iterated_local_search(&problem, 20, iterations, rng_seed);
        rng_seed += 1;

        let (_, _, positions_iterated_30, _) = heuristic_iterated_local_search(&problem, 30, iterations, rng_seed);
        rng_seed += 1;

        let (_, _, positions_iterated_40, _) = heuristic_iterated_local_search(&problem, 40, iterations, rng_seed);
        rng_seed += 1;

        let (_, _, positions_iterated_50, _) = heuristic_iterated_local_search(&problem, 50, iterations, rng_seed);
        rng_seed += 1;

        let (_, _, positions_iterated_75, _) = heuristic_iterated_local_search(&problem, 75, iterations, rng_seed);
        rng_seed += 1;

        println!("{} || {:<12.5e} || {:<12.5e} | {:<12.5e} | {:<12.5e} || {:<12.5e} | {:<12.5e} | {:<12.5e} || {:<12.5e} | {:<12.5e} | {:<12.5e} | {:<12.5e} | {:<12.5e} | {:<12.5e} | {:<12.5e}",
                 file,
                 // Exact
                 compute_error(&problem, &positions_exact),
                 // Greedy
                 compute_error(&problem, &positions_greedy_one),
                 compute_error(&problem, &positions_greedy_two),
                 compute_error(&problem, &positions_greedy_worst),
                 // Local search
                 compute_error(&problem, &positions_local_forward),
                 compute_error(&problem, &positions_local_backward),
                 compute_error(&problem, &positions_local_best),
                 // Iterated local search
                 compute_error(&problem, &positions_iterated_5),
                 compute_error(&problem, &positions_iterated_10),
                 compute_error(&problem, &positions_iterated_20),
                 compute_error(&problem, &positions_iterated_30),
                 compute_error(&problem, &positions_iterated_40),
                 compute_error(&problem, &positions_iterated_50),
                 compute_error(&problem, &positions_iterated_75),
        );
    }

    println!("\n1ptq || 5%           | 10%          | 20%          | 30%          | 40%          | 50%           ");

    let problem = load_problem("1ptq");
    let ils_seed = 0;

    let (_, _, _, errors_iterated_5) = heuristic_iterated_local_search2(&problem, 5, 1000, ils_seed);
    let (_, _, _, errors_iterated_10) = heuristic_iterated_local_search2(&problem, 10, 1000, ils_seed + 1);
    let (_, _, _, errors_iterated_20) = heuristic_iterated_local_search2(&problem, 20, 1000, ils_seed + 2);
    let (_, _, _, errors_iterated_30) = heuristic_iterated_local_search2(&problem, 30, 1000, ils_seed + 3);
    let (_, _, _, errors_iterated_40) = heuristic_iterated_local_search2(&problem, 40, 1000, ils_seed + 4);
    let (_, _, _, errors_iterated_50) = heuristic_iterated_local_search2(&problem, 50, 1000, ils_seed + 5);
    let (_, _, _, errors_iterated_75) = heuristic_iterated_local_search2(&problem, 75, 1000, ils_seed + 6);

    for i in 0..1000 {
        println!("{:<5} || {:<12.5e} | {:<12.5e} | {:<12.5e} | {:<12.5e} | {:<12.5e} | {:<12.5e} | {:<12.5e}",
                 i,
                 errors_iterated_5[i],
                 errors_iterated_10[i],
                 errors_iterated_20[i],
                 errors_iterated_30[i],
                 errors_iterated_40[i],
                 errors_iterated_50[i],
                 errors_iterated_75[i],
        );
    }
}
