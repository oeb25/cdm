use itertools::Itertools;
use tracing::debug;

use crate::{Natural, Polynomial, Ring};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PrimitiveRootOfUnity<R> {
    nth: Natural,
    value: R,
    inverse: R,
}

impl<R> PrimitiveRootOfUnity<R> {
    /// Definition 8.5
    pub fn new(n: Natural, omega: R) -> Option<Self>
    where
        R: Ring,
    {
        let n = n.into();
        let mut inverse = omega.clone();
        let mut pow = omega.clone();

        for i in 1..u128::from(n) {
            // println!("{omega:?}^{i} = {pow:?}");
            if pow.is_one() || pow.is_zero() {
                return None;
            }
            inverse = pow.clone();
            pow = pow * omega.clone();
        }
        // println!("{omega:?}^{n:?} = {pow:?}");
        if !pow.is_one() {
            return None;
        }

        Some(Self {
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
    pub fn sq(self) -> Self
    where
        R: Clone + std::ops::Mul<Output = R>,
    {
        Self {
            nth: self.nth * 2,
            value: self.value.clone() * self.value,
            inverse: self.inverse.clone() * self.inverse,
        }
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

/// Algorithm 8.14 Fast Fourier Transform (FFT)
pub fn fft<R: Ring + std::fmt::Debug>(
    k: Natural,
    f: Polynomial<R>,
    omega: PrimitiveRootOfUnity<R>,
) -> Vec<R> {
    // pub fn fft<R: Ring + std::fmt::Debug>(k: Natural, f: impl DFTInput<R>, omega: R) -> Vec<R> {

    let span = tracing::span!(
        tracing::Level::DEBUG,
        "FFT",
        omega = format!("{:?}", omega.clone().inner()),
        f = format!("{:?}", f)
    );
    let _enter = span.enter();

    if k == 0 {
        return vec![f.evaluate_at(R::zero())];
    }
    let n_div_2 = 2u128.pow((k - 1) as _);

    let r0 = f
        .iter()
        .take(n_div_2 as _)
        .map(|(c, pow)| c.clone() + f.coef_at(n_div_2 + pow))
        .collect_vec();
    let r0 = Polynomial::new(r0);
    let r1 = f
        .iter()
        .take(n_div_2 as _)
        .map(|(c, pow)| (c.clone() - f.coef_at(n_div_2 + pow)) * omega.pow(pow))
        .collect_vec();
    let r1 = Polynomial::new(r1);

    debug!(
        "r0 = {:?} = {r0:?}",
        f.iter()
            .take(n_div_2 as _)
            .map(|(c, pow)| format!("({c:?} + {:?})*x^{pow}", f.coef_at(n_div_2 + pow)))
            .collect_vec()
    );
    debug!(
        "r1 = {:?} = {r1:?}",
        f.iter()
            .take(n_div_2 as _)
            .map(|(c, pow)| format!(
                "({c:?} - {:?})*{:?}*x^{pow}",
                f.coef_at(n_div_2 + pow),
                omega.pow(pow)
            ))
            .collect_vec()
    );

    let omega_sq =
        PrimitiveRootOfUnity::new(omega.n() / 2, omega.clone().inner() * omega.inner()).unwrap();

    let r0_eval = fft(k - 1, r0, omega_sq.clone());
    let r1_eval = fft(k - 1, r1, omega_sq);

    let res = r0_eval
        .into_iter()
        .interleave(r1_eval.into_iter())
        .collect_vec();

    assert_eq!(res.len(), 2usize.pow(u128::from(k) as u32));

    debug!("result = {res:?}");

    res
}

#[cfg(test)]
mod tests {
    use crate::{dft::PrimitiveRootOfUnity, Finite, Polynomial};

    #[test]
    fn simple_fft() {
        type R = Finite<5>;

        let f = Polynomial::new([1, 1, 0, 1u128].map(R::from).to_vec());

        eprintln!("{f:?}");

        let res = super::fft(2, f, PrimitiveRootOfUnity::new(4, 2u128.into()).unwrap());

        // panic!("{res:?}");
    }
}
