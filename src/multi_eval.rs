use itertools::Itertools;

use crate::{field::Field, Group, Polynomial, Ring};

/// ALGORITHM 10.5 Going down the subproduct tree.
pub fn multi_eval<F: Ring + std::fmt::Debug>(p: Polynomial<F>, us: &[F]) -> Vec<F> {
    eprintln!("Evaluate {p:?} at {us:?}");

    let n = us.len();
    let k = n.ilog2() as usize;

    if n == 1 {
        return vec![p.evaluate_at(F::zero())];
    }

    let m = building_up_the_subproduct_tree(us);

    let r0 = p.clone() % m[k - 1][0].clone();
    let r1 = p % m[k - 1][1].clone();

    eprintln!("r0 = {r0:?}");
    eprintln!("r1 = {r1:?}");

    multi_eval(r0, &us[0..n / 2])
        .into_iter()
        .chain(multi_eval(r1, &us[n / 2..]))
        .collect()
}

pub use linear_combination_for_linear_moduli as algorithm_10_9;
/// ALGORITHM 10.9 Linear combination for linear moduli.
pub fn linear_combination_for_linear_moduli<F: Field + std::fmt::Debug>(
    us: &[F],
    cs: &[F],
) -> Polynomial<F> {
    let n = us.len();
    let k = n.ilog2() as usize;

    if n == 1 {
        return Polynomial::new(vec![cs[0].clone()]);
    }

    let s0 = linear_combination_for_linear_moduli(&us[0..n / 2], &cs[0..n / 2]);
    let s1 = linear_combination_for_linear_moduli(&us[n / 2..], &cs[n / 2..]);

    let m = building_up_the_subproduct_tree(us);

    &m[k - 1][1] * &s0 + &m[k - 1][0] * &s1
}

/// ALGORITHM 10.3 Building up the subproduct tree.
pub fn building_up_the_subproduct_tree<F: Ring + std::fmt::Debug>(
    us: &[F],
) -> Vec<Vec<Polynomial<F>>> {
    let n = us.len();
    let k = n.ilog2() as usize;

    let mut m: Vec<Vec<Polynomial<F>>> = vec![vec![]];
    for u in us {
        m[0].push(Polynomial::x() - u);
    }
    for i in 1..=k {
        m.push(vec![
            Polynomial::new(vec![F::zero()]);
            2usize.pow((k - i) as _)
        ]);
        for j in 0..=2usize.pow((k - i) as _) - 1 {
            m[i][j] = &m[i - 1][2 * j] * &m[i - 1][2 * j + 1];
        }
    }
    m
}

/// ALGORITHM 10.11 Fast interpolation.
pub fn fast_interpolation<F: Field + std::fmt::Debug>(us: &[F], vs: &[F]) -> Polynomial<F> {
    let n = us.len();
    let k = n.ilog2() as usize;

    // 1.
    let matrix = building_up_the_subproduct_tree(us);

    // 2.
    let m = matrix[k][0].clone();
    eprintln!("m = {m:?}");
    let f = m.diff();
    eprintln!("m' = {f:?}");
    let m_diff_eval = multi_eval(f, us);
    let s = m_diff_eval
        .iter()
        .zip(vs)
        .map(|(v_diff, v)| v.clone() / v_diff.clone())
        .collect_vec();

    // 3.
    linear_combination_for_linear_moduli(us, &s)
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;
    use crate::{Polynomial, Rational};

    #[test]
    fn multi_example() {
        let f = Polynomial::new([1, 2, 3, 4, 5, 6, 7, 8].map(Rational::from).to_vec());
        let us = (-3..=4).map(Rational::from).collect_vec();

        let eval = multi_eval(f, &us);

        panic!("{eval:#?}");
    }

    #[test]
    fn interpolation() {
        let us = [0, 1, 2, 3].map(Rational::from);
        let vs = [1, 2, 4, 8].map(Rational::from);
        let result = fast_interpolation(&us, &vs);

        panic!("{result:?}");
    }
}
