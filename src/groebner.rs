use itertools::Itertools;

use crate::{
    ch21::multivariate_division_with_remainder, field::Field, mono::MonomialOrder,
    multivariate_polynomials::MultivariatePolynomial, Group, Ring,
};

#[derive(Debug)]
pub struct Basis<F: Ring, O: MonomialOrder<F>> {
    gs: Vec<MultivariatePolynomial<F, O>>,
}

pub fn buchbergers_algorithm<F, O>(fs: &[MultivariatePolynomial<F, O>]) -> Basis<F, O>
where
    F: Field,
    O: MonomialOrder<F>,
{
    let mut basis = fs.to_vec();
    basis.dedup();

    loop {
        let mut new = vec![];

        for (i, f1) in basis.iter().enumerate() {
            for f2 in &basis[i + 1..] {
                let s = f1.s_polynomial(f2);
                let result = multivariate_division_with_remainder(&s, &basis);

                let r = result.remainder();

                if !r.is_zero() && !new.contains(&r) && !basis.contains(&r) {
                    if !multivariate_division_with_remainder(&r, &new)
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

pub fn minimize_groebner_basis<F, O>(basis: &mut Basis<F, O>)
where
    F: Field,
    O: MonomialOrder<F>,
{
    let mut i = 0;

    while i < basis.gs.len() {
        basis.gs.swap(0, i);

        let f = basis.gs[0].leading_term().into();
        let rest = basis.gs[1..]
            .iter()
            .map(|g| g.leading_term().into())
            .collect_vec();

        // let r = multivariate_division_with_remainder(ord, &basis.gs[0], &basis.gs[1..]).remainder();
        let r = multivariate_division_with_remainder(&f, &rest).remainder();

        if r.is_zero() {
            basis.gs.remove(0);
        } else {
            i += 1;
        }
    }

    for g in &mut basis.gs {
        *g = g.minimize();
    }
}
