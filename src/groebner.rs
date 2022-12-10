use itertools::Itertools;

use crate::{
    field::Field,
    identity::{Addition, Identity, Multiplication},
    mono::MonomialOrder,
    multivariate_polynomials::{multivariate_division_with_remainder, MultivariatePolynomial},
    Group,
};

#[derive(Debug)]
pub struct Basis<F: Identity<Addition> + Identity<Multiplication>> {
    gs: Vec<MultivariatePolynomial<F>>,
}

pub fn buchbergers_algorithm<F>(
    ord: &impl MonomialOrder<F>,
    fs: &[MultivariatePolynomial<F>],
) -> Basis<F>
where
    F: Field,
{
    let mut basis = fs.to_vec();
    basis.dedup();

    loop {
        let mut new = vec![];

        for (i, f1) in basis.iter().enumerate() {
            for f2 in &basis[i + 1..] {
                let s = f1.s_polynomial(f2, ord);
                let result = multivariate_division_with_remainder(ord, &s, &basis);

                let mut r = result.remainder();

                if !r.is_zero() && !new.contains(&r) && !basis.contains(&r) {
                    r.sort_by(ord);

                    if !multivariate_division_with_remainder(ord, &r, &new)
                        .remainder()
                        .is_zero()
                    {
                        new.push(r);
                    }
                }
            }
        }

        if new.is_empty() {
            break;
        }

        basis.extend_from_slice(&new);
    }

    Basis { gs: basis }
}

pub fn minimize_groebner_basis<F>(ord: &impl MonomialOrder<F>, basis: &mut Basis<F>)
where
    F: Field,
{
    let mut i = 0;

    while i < basis.gs.len() {
        basis.gs.swap(0, i);

        let f = basis.gs[0].leading_term(ord).into();
        let rest = basis.gs[1..]
            .iter()
            .map(|g| g.leading_term(ord).into())
            .collect_vec();

        // let r = multivariate_division_with_remainder(ord, &basis.gs[0], &basis.gs[1..]).remainder();
        let r = multivariate_division_with_remainder(ord, &f, &rest).remainder();

        if r.is_zero() {
            basis.gs.remove(0);
        } else {
            i += 1;
        }
    }

    for g in &mut basis.gs {
        *g = g.normalized().minimize(ord);
    }
}
