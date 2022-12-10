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
    assert_eq!(2u128.pow(k as _), omega.n());

    let alpha = dft::fft(k, f, omega.clone());
    let beta = dft::fft(k, g, omega.clone());
    eprintln!("alpha = {alpha:?}");
    eprintln!("beta = {beta:?}");
    let gamma = alpha
        .into_iter()
        .zip(beta)
        .map(|(l, r)| l * r)
        .collect_vec();
    eprintln!("gamma = {gamma:?}");

    let n_in_r = (0..omega.n() as _)
        .map(|_| R::one())
        .fold(R::zero(), |a, b| a + b);

    Polynomial::new(dft::fft(
        k,
        Polynomial::new(gamma),
        PrimitiveRootOfUnity::new(omega.n(), omega.inverse()).unwrap(),
    )) * Polynomial::new(vec![n_in_r.multiplicative_inverse().unwrap()])
}
