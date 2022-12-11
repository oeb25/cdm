pub mod ch08;
pub mod count_ops;
pub mod cra;
pub mod dft;
pub mod domain;
pub mod euclidean_domain;
pub mod field;
pub mod finite;
pub mod gaussian_integers;
pub mod groebner;
pub mod group;
pub mod identity;
pub mod integers;
pub mod latex;
pub mod mono;
pub mod multi_eval;
pub mod multivariate_polynomials;
pub mod naturals;
pub mod newton_interpolation;
pub mod polynomials;
pub mod rationals;
pub mod reals;
pub mod ring;

pub use finite::Finite;
pub use group::Group;
pub use integers::Integer;
pub use mono::Monomial;
pub use naturals::Natural;
pub use polynomials::Polynomial;
pub use rationals::Rational;
pub use reals::Real;
pub use ring::Ring;

use itertools::Itertools;

use crate::{
    dft::PrimitiveRootOfUnity,
    euclidean_domain::ExtendedEuclideanAlgorithm,
    gaussian_integers::Gaussian,
    latex::ToLatex,
    mono::{MonomialOrder, PLex},
    multivariate_polynomials::{multivariate_division_with_remainder, MultivariatePolynomial},
};

pub const fn digits_superscript(c: char) -> Option<char> {
    Some(match c {
        '0' => '⁰',
        '1' => '¹',
        '2' => '²',
        '3' => '³',
        '4' => '⁴',
        '5' => '⁵',
        '6' => '⁶',
        '7' => '⁷',
        '8' => '⁸',
        '9' => '⁹',
        '+' => '⁺',
        '-' => '⁻',
        '=' => '⁼',
        '(' => '⁽',
        ')' => '⁾',
        _ => return None,
    })
}
pub const fn digits_subscript(c: char) -> Option<char> {
    Some(match c {
        '0' => '₀',
        '1' => '₁',
        '2' => '₂',
        '3' => '₃',
        '4' => '₄',
        '5' => '₅',
        '6' => '₆',
        '7' => '₇',
        '8' => '₈',
        '9' => '₉',
        '+' => '₊',
        '-' => '₋',
        '=' => '₌',
        '(' => '₍',
        ')' => '₎',
        _ => return None,
    })
}
pub fn num_to_superscript(n: i128) -> String {
    n.to_string().chars().flat_map(digits_superscript).collect()
}
pub fn num_to_subscript(n: i128) -> String {
    n.to_string().chars().flat_map(digits_subscript).collect()
}

pub fn mpoly<F: Ring, const N: usize>(
    terms: [(impl Into<F>, &[u128]); N],
) -> MultivariatePolynomial<F> {
    let terms = terms
        .into_iter()
        .map(|(c, ps)| Monomial::new(c.into(), ps.to_vec()))
        .collect();

    MultivariatePolynomial { terms }
}
pub fn mpoly_rat<const N: usize>(terms: [(f64, &[u128]); N]) -> MultivariatePolynomial<Rational> {
    let terms = terms
        .into_iter()
        .map(|(c, ps)| Monomial::new(Rational::approximate(c), ps.to_vec()))
        .collect();

    MultivariatePolynomial { terms }
}
