//! # Newton iteration

use tracing::debug;

use crate::{Group, Natural, Polynomial, Ring};

/// Algorithm 9.3 Inversion using Newton iteration
///
/// - Input: f ∈ D[x] with f(0) = 1, and l ∈ N.
/// - Output: g ∈ D[x] satisfying fg ≡ 1 mod x^l.
pub fn inversion_newton_iteration<D: Ring>(f: Polynomial<D>, l: Natural) -> Polynomial<D> {
    let r = (l as f64).log2().ceil() as Natural;

    let scope = tracing::span!(
        tracing::Level::DEBUG,
        "Inversion Newton iteration",
        f = format!("{f:?}"),
        l = l.to_string()
    );
    let _enter = scope.enter();

    let mut g: Polynomial<D> = Polynomial::new(vec![D::one()]);
    for i in 1..=r {
        g = g.clone() + g.clone()
            - f.clone() * (g.clone() * g.clone()) % Polynomial::one().times_x(2.pow(i) as _);
        debug!("g_{i} = {g:?}");
    }
    assert_eq!((f * g.clone()) % Polynomial::one().times_x(l), Ring::one());
    g
}

