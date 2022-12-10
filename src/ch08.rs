use itertools::Itertools;

use crate::{
    dft::{self, PrimitiveRootOfUnity},
    Natural, Polynomial, Ring,
};

/// Algorithm 8.16 Fast convolution
pub fn fast_convolution<R: Ring>(
    k: Natural,
    f: Polynomial<R>,
    g: Polynomial<R>,
    omega: PrimitiveRootOfUnity<R>,
) -> Polynomial<R> {
    assert_eq!(Natural::from(2).pow(k), omega.n());

    let alpha = dft::fft(k, f, omega.clone().inner());
    let beta = dft::fft(k, g, omega.clone().inner());
    eprintln!("alpha = {alpha:?}");
    eprintln!("beta = {beta:?}");
    let gamma = alpha
        .into_iter()
        .zip(beta)
        .map(|(l, r)| l * r)
        .collect_vec();
    eprintln!("gamma = {gamma:?}");

    let n_in_r = (0..i128::from(omega.n()))
        .map(|_| R::one())
        .fold(R::zero(), |a, b| a + b);

    Polynomial::new(dft::fft(k, Polynomial::new(gamma), omega.inverse()))
        * Polynomial::new(vec![n_in_r.multiplicative_inverse().unwrap()])
}
