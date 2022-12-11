//! # Newton iteration

use tracing::debug;

use crate::{Group, Natural, Polynomial, Ring};

/// Algorithm 9.3 Inversion using Newton iteration.
///
/// - Input: `f ∈ D[x]` with `f(0) = 1`, and `l ∈ N`.
/// - Output: `g ∈ D[x]` satisfying `fg ≡ 1 mod x^l`.
pub fn inversion_newton_iteration<D: Ring>(f: Polynomial<D>, l: Natural) -> Polynomial<D> {
    let r = (l as f64).log2().ceil() as Natural;

    let scope = tracing::span!(
        tracing::Level::DEBUG,
        "(9.3) Inversion Newton iteration",
        f = format!("{f:?}"),
        l = l.to_string()
    );
    let _enter = scope.enter();

    let mut g: Polynomial<D> = Polynomial::new(vec![D::one()]);
    for i in 1..=r {
        g = (g.clone() + g.clone() - f.clone() * (g.clone() * g.clone())).rem_pow(2.pow(i) as _);
        debug!("g_{i} = {g:?}");
    }
    assert_eq!((f * g.clone()) % Polynomial::one().times_x(l), Ring::one());
    debug!("g = {g:?}");
    g
}

/// Algorithm 9.5 Fast division with remainder.
pub fn fast_division_with_remainder<D: Ring>(
    a: Polynomial<D>,
    b: Polynomial<D>,
) -> (Polynomial<D>, Polynomial<D>) {
    let scope = tracing::span!(
        tracing::Level::DEBUG,
        "(9.5) Fast division with remainder",
        a = format!("{a:?}"),
        b = format!("{b:?}"),
    );
    let _enter = scope.enter();

    assert!(b != Polynomial::zero() && b.is_monic());

    if a.deg() < b.deg() {
        return (Polynomial::zero(), a);
    }
    let m = a.deg() - b.deg();

    let q_star =
        (a.rev(a.deg()) * inversion_newton_iteration(b.rev(b.deg()), m + 1)).rem_pow(m + 1);
    debug!("q_star = {q_star:?}");

    let q = q_star.rev(m);
    let r = a.clone() - b.clone() * q.clone();
    debug!("q = {q:?}, r = {r:?}");
    (q, r)
}

/// Algorithm 9.10 p-adic inversion using Newton iteration.
pub fn p_adic_inversion_using_newton_iteration() {
    unimplemented!()
}
