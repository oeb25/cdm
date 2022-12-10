use itertools::Itertools;

use crate::{Natural, Polynomial, Ring};


#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PrimitiveRootOfUnity<R> {
    nth: Natural,
    value: R,
    inverse: R,
}

impl<R> PrimitiveRootOfUnity<R> {
    /// Definition 8.5
    pub fn new(n: impl Into<Natural>, omega: R) -> Option<PrimitiveRootOfUnity<R>>
    where
        R: Ring,
    {
        let n = n.into();
        let mut inverse = omega.clone();
        let mut pow = omega.clone();

        for _ in 1..u128::from(n) {
            if pow.is_one() || pow.is_zero() {
                return None;
            }
            inverse = pow.clone();
            pow = pow * omega.clone();
        }
        if !pow.is_one() {
            return None;
        }

        Some(PrimitiveRootOfUnity {
            nth: n,
            value: omega,
            inverse,
        })
    }
    pub fn n(&self) -> Natural {
        self.nth
    }
    pub fn inner(self) -> R {
        self.value
    }
    pub fn inverse(self) -> R {
        self.inverse
    }
}
impl<R> std::ops::Deref for PrimitiveRootOfUnity<R> {
    type Target = R;

    fn deref(&self) -> &Self::Target {
        &self.value
    }
}


pub trait DFTInput<R> {
    fn coef_at(&self, pow: Natural) -> R;
    fn evaluate_at(&self, x: R) -> R;
    fn elements(&self) -> Vec<(&R, Natural)>;
}

impl<R> DFTInput<R> for Polynomial<R>
where
    R: Ring,
{
    fn coef_at(&self, pow: Natural) -> R {
        self.coef_at(pow)
    }

    fn evaluate_at(&self, x: R) -> R {
        self.evaluate_at(x)
    }

    fn elements(&self) -> Vec<(&R, Natural)> {
        self.iter().collect()
    }
}

pub fn fft<R: Ring + std::fmt::Debug>(k: Natural, f: Polynomial<R>, omega: R) -> Vec<R> {
    // pub fn fft<R: Ring + std::fmt::Debug>(k: Natural, f: impl DFTInput<R>, omega: R) -> Vec<R> {
    if k.is_zero() {
        vec![f.evaluate_at(R::zero())]
    } else {
        let n_div_2 = Natural::from(2).pow(k - 1.into());

        let r0 = f
            .iter()
            .take(n_div_2.into())
            .map(|(c, pow)| c.clone() + f.coef_at(n_div_2 + pow))
            .collect_vec();
        let r0 = Polynomial::new(r0);
        let r1 = f
            .iter()
            .take(n_div_2.into())
            .map(|(c, pow)| (c.clone() - f.coef_at(n_div_2 + pow)) * omega.pow(pow))
            .collect_vec();
        let r1 = Polynomial::new(r1);

        // eprintln!("{r0:?}");
        // eprintln!("{r1:?}");

        let omega_sq = omega.clone() * omega;

        let r0_eval = fft(k - 1.into(), r0, omega_sq.clone());
        let r1_eval = fft(k - 1.into(), r1, omega_sq);

        let res = r0_eval
            .into_iter()
            .interleave(r1_eval.into_iter())
            .collect_vec();

        assert_eq!(res.len(), 2usize.pow(u128::from(k) as u32));

        res
    }
}

#[cfg(test)]
mod tests {
    use crate::{Finite, Polynomial};

    #[test]
    fn simple_fft() {
        type R = Finite<5>;

        let f: Polynomial<R> = Polynomial::new([1, 1, 0, 1u128].map(R::from).to_vec());

        eprintln!("{f:?}");

        let res = super::fft(2.into(), f, 2u128.into());

        // panic!("{res:?}");
    }
}
