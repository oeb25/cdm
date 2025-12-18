//! # Fast multiplication

use itertools::Itertools;
use tracing::debug;

use crate::{
    dft::{self, PrimitiveRootOfUnity},
    field::Field,
    Natural, Polynomial, Ring,
};

/// Algorithm 8.1 Karatsuba’s polynomial multiplication algorithm
pub fn karatsubas_polynomial_multiplication_algorithm<R: Field>(
    k: Natural,
    f: Polynomial<R>,
    g: Polynomial<R>,
) -> Polynomial<R> {
    let scope = tracing::span!(
        tracing::Level::DEBUG,
        "(8.1) Karatsuba",
        k = format!("{k:?}"),
        f = format!("{f:?}"),
        g = format!("{g:?}"),
    );
    let _enter = scope.enter();

    let n = 2u128.pow(k as _);
    assert!(f.deg() < n && g.deg() < n);
    if n == 1 {
        return f * g;
    }
    let f1 = Polynomial::new(
        f.iter()
            .skip(n as usize / 2)
            .map(|(c, _)| c.clone())
            .collect(),
    );
    let f0 = Polynomial::new(
        f.iter()
            .take(n as usize / 2)
            .map(|(c, _)| c.clone())
            .collect(),
    );

    let g1 = Polynomial::new(
        g.iter()
            .skip(n as usize / 2)
            .map(|(c, _)| c.clone())
            .collect(),
    );
    let g0 = Polynomial::new(
        g.iter()
            .take(n as usize / 2)
            .map(|(c, _)| c.clone())
            .collect(),
    );

    // debug!("f0 = {f0:?}");
    // debug!("f1 = {f1:?}");
    // debug!("g0 = {g0:?}");
    // debug!("g1 = {g1:?}");

    let fg0 = karatsubas_polynomial_multiplication_algorithm(k - 1, f0.clone(), g0.clone());
    let fg1 = karatsubas_polynomial_multiplication_algorithm(k - 1, f1.clone(), g1.clone());
    let mixed = karatsubas_polynomial_multiplication_algorithm(k - 1, f0 + f1, g0 + g1);

    debug!("fg0 = {fg0:?}");
    debug!("fg1 = {fg1:?}");
    debug!("mixed = {mixed:?}");

    let res = fg1.times_x(n) + (mixed - fg0.clone() - fg1).times_x(n / 2) + fg0;

    debug!("result = {res:?}");

    res
}

/// Algorithm 8.14 Fast Fourier Transform (FFT)
pub use dft::fft as fast_fourier_transform;

/// Algorithm 8.16 Fast convolution
pub fn fast_convolution<R: Ring>(
    k: Natural,
    f: Polynomial<R>,
    g: Polynomial<R>,
    omega: PrimitiveRootOfUnity<R>,
) -> Polynomial<R> {
    let scope = tracing::span!(
        tracing::Level::DEBUG,
        "(8.16) Fast convolution",
        f = format!("{f:?}"),
        g = format!("{g:?}"),
        omega = format!("{:?}", omega.clone().inner()),
    );
    let _enter = scope.enter();

    assert_eq!(2u128.pow(k as _), omega.n());

    let alpha = dft::fft(k, f, omega.clone());
    let beta = dft::fft(k, g, omega.clone());
    debug!("alpha = {alpha:?}");
    debug!("beta = {beta:?}");
    let gamma = alpha
        .into_iter()
        .zip(beta)
        .map(|(l, r)| l * r)
        .collect_vec();
    debug!("gamma = {gamma:?}");

    let n_in_r = (0..omega.n() as _)
        .map(|_| R::one())
        .fold(R::zero(), |a, b| a + b);

    Polynomial::new(dft::fft(
        k,
        Polynomial::new(gamma),
        PrimitiveRootOfUnity::new(omega.n(), omega.inverse()).unwrap(),
    )) * Polynomial::new(vec![n_in_r.multiplicative_inverse().unwrap()])
}

/// Algorithm 8.20 Fast negative wrapped convolution
pub fn fast_negative_wrapped_convolution() {
    unimplemented!()
}
