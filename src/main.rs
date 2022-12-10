mod ch08;
pub mod count_ops;
pub mod dft;
mod domain;
mod euclidean_domain;
mod field;
pub mod finite;
mod gaussian_integers;
mod groebner;
mod group;
mod identity;
mod integers;
mod latex;
mod mono;
mod multi_eval;
mod multivariate_polynomials;
mod naturals;
mod newton_interpolation;
mod polynomials;
mod rationals;
mod reals;
mod ring;

pub use finite::Finite;
pub use group::Group;
pub use integers::Integer;
use mono::Monomial;
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

fn main() {
    // test_exam();

    type R = Integer;

    let res = ExtendedEuclideanAlgorithm::perform(&R::from(12i128), &R::from(29i128));

    println!("{res:?}")
}

fn exercise_8_10() {
    type R = Finite<17>;

    let f = Polynomial::new([3, -4, 3, 5i128].map(R::from).to_vec());
    let g = Polynomial::new([-2, 7, -5, 2i128].map(R::from).to_vec());

    // (i)
    let omega = PrimitiveRootOfUnity::new(8, R::from(2i128)).expect("it is");

    // (ii)
    let h = f.clone() * g.clone();

    // (iii)
    for j in 0..8 {
        let o = omega.pow(Natural::from(j));

        let alpha_j = f.evaluate_at(o);
        let beta_j = g.evaluate_at(o);

        let gamma_j = alpha_j * beta_j;
        dbg!(gamma_j);
        dbg!(h.evaluate_at(o));
    }

    // (v)
    let res = ch08::fast_convolution(3.into(), f, g, omega);
    println!("{:?}", h);
    println!("{:?}", res);
}

fn test_exam() {
    fn question_1() {
        let f = Gaussian::new(7, 8).map(Integer::from);
        let g = Gaussian::new(2, 3).map(Integer::from);

        let result = ExtendedEuclideanAlgorithm::perform(&f, &g);
        println!("{result:#?}");

        dbg!(result.gcd());
    }

    question_1();
}

fn old_exam() {
    struct FuckOrdering;

    impl MonomialOrder<Finite<5>> for FuckOrdering {
        fn ord(&self, l: &Monomial<Finite<5>>, r: &Monomial<Finite<5>>) -> std::cmp::Ordering {
            let a = l
                .powers()
                .iter()
                .enumerate()
                .map(|(pow, n)| Natural::from(pow as u128 + 1) * *n)
                .fold(Natural::zero(), |a, b| a + b);
            let b = r
                .powers()
                .iter()
                .enumerate()
                .map(|(pow, n)| Natural::from(pow as u128 + 1) * *n)
                .fold(Natural::zero(), |a, b| a + b);

            b.cmp(&a).then_with(|| PLex::default().ord(l, r))
        }
    }

    type R = Finite<5>;

    let ord = FuckOrdering;

    let g1: MultivariatePolynomial<R> = mpoly([(1i128, &[3]), (-1i128, &[0, 0, 1])]);
    let g2: MultivariatePolynomial<R> = mpoly([(-1i128, &[1]), (1i128, &[0, 1])]);

    eprintln!("{:?}", g1.leading_term(&ord));
    eprintln!("{:?}", g2.leading_term(&ord));

    let res = multivariate_division_with_remainder(
        &ord,
        &g1.s_polynomial(&g2, &ord),
        &[g1.clone(), g2.clone()],
    );

    eprintln!("{}", res.ascii_table(&ord));

    let mut G = groebner::buchbergers_algorithm(&ord, &[g1, g2]);

    eprintln!("{G:?}");
    groebner::minimize_groebner_basis(&ord, &mut G);
    eprintln!("{G:?}");
}

fn home_work2() {
    // let f: Integer = 126.into();
    // let g: Integer = 35.into();

    // let res = ExtendedEuclideanAlgorithm::perform(&f, &g);
    // dbg!(&res);
    // dbg!(res.gcd());

    // let f =
    //     Polynomial::<Rational>::new(vec![18, -42, 30, -6].into_iter().rev().map_into().collect());
    // let g = Polynomial::new(vec![-12, 10, -2].into_iter().rev().map_into().collect());
    // let res = ExtendedEuclideanAlgorithm::perform(&f, &g);
    // println!("{:#?}", res);
    // println!("{:?}", res.gcd());

    let ord = PLex(vec![2, 0, 1]);

    let f1: MultivariatePolynomial<Rational> = MultivariatePolynomial {
        terms: vec![
            Monomial::new(1, vec![2]),
            Monomial::new(1, vec![0, 2]),
            Monomial::constant(-1),
        ],
    }
    .normalized()
    .sorted_by(&ord);
    let f2: MultivariatePolynomial<Rational> = MultivariatePolynomial {
        terms: vec![
            Monomial::new(2, vec![1, 1]),
            Monomial::new(-2, vec![1, 0, 1]),
            Monomial::new(-2, vec![0, 1]),
        ],
    }
    .normalized()
    .sorted_by(&ord);
    let f3: MultivariatePolynomial<Rational> = MultivariatePolynomial {
        terms: vec![
            Monomial::new(1, vec![2]),
            Monomial::new(-2, vec![0, 1, 1]),
            Monomial::new(-2, vec![1]),
            Monomial::constant(1),
        ],
    }
    .normalized()
    .sorted_by(&ord);

    eprintln!("{f1:?}");
    eprintln!("{f2:?}");
    // eprintln!("{f3:?}");

    let s12 = f1.s_polynomial(&f2, &ord);

    // -2*x^2*y - 2*y^2*z + 2*x*y + 2*z
    eprintln!("s12 = {s12:?}");
    let s12 = s12
        * MultivariatePolynomial {
            terms: vec![Monomial::constant(1)],
        };
    eprintln!("s12 = {s12:?}");

    let fs = [f1, f2, f3];

    let table = multivariate_division_with_remainder::<Rational>(&ord, &s12, &fs);

    // println!("{}", table.latex_table(&ord));

    let mut res = groebner::buchbergers_algorithm(&ord, &fs);

    eprintln!("One: {res:#?}");
    groebner::minimize_groebner_basis(&ord, &mut res);
    eprintln!("Two: {res:#?}");

    let g1 = mpoly_rat([(1., &[0, 0, 1]), (-9., &[0, 5, 0]), (1. / 2., &[0, 3, 0])]);
    let g2 = mpoly_rat([
        (1., &[1, 0, 0]),
        (9. / 2., &[0, 4, 0]),
        (1. / 2., &[0, 2, 0]),
        (-1., &[]),
    ]);
    let g3 = mpoly_rat([(1., &[0, 6, 0]), (-5. / 9., &[0, 4])]);

    eprintln!("{g1:?}");
    eprintln!("{g2:?}");
    eprintln!("{g3:?}");

    let gs = [g1, g2, g3];
    let mut basis = groebner::buchbergers_algorithm(&ord, &gs);

    eprintln!("Pre:  {basis:#?}");
    groebner::minimize_groebner_basis(&ord, &mut basis);
    eprintln!("Post: {basis:#?}");

    for i in 0..gs.len() {
        for j in i + 1..gs.len() {
            let s = gs[i].s_polynomial(&gs[j], &ord);
            let r = multivariate_division_with_remainder(&ord, &s, &gs);

            println!("{}", r.latex_table(&ord));
        }
    }

    eprintln!("{}", gs.iter().map(|g| g.to_latex()).format("\n"));
}

fn mpoly<F: Ring, const N: usize>(
    terms: [(impl Into<F>, &[u128]); N],
) -> MultivariatePolynomial<F> {
    let terms = terms
        .into_iter()
        .map(|(c, ps)| Monomial::new(c.into(), ps.to_vec()))
        .collect();

    MultivariatePolynomial { terms }
}
fn mpoly_rat<const N: usize>(terms: [(f64, &[u128]); N]) -> MultivariatePolynomial<Rational> {
    let terms = terms
        .into_iter()
        .map(|(c, ps)| Monomial::new(Rational::approximate(c), ps.to_vec()))
        .collect();

    MultivariatePolynomial { terms }
}
