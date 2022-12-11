use cdm::{
    ch21::{buchbergers_algorithm, minimize_groebner_basis, multivariate_division_with_remainder},
    dft::{fft, PrimitiveRootOfUnity},
    mono::{plex, MonomialOrder},
    multivariate_polynomials::MultivariatePolynomial,
    Finite, Monomial, Natural, Polynomial,
};
use tracing::{debug, info};

fn main() {
    tracing_subscriber::fmt::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .without_time()
        .init();

    q1();
    q2();
    q3();
}

#[tracing::instrument]
fn q1() {
    info!("Question 1");

    type R = Finite<5>;
    let n = 4u128;

    // (a)
    let omega = PrimitiveRootOfUnity::new(n, R::from(2i128)).unwrap();
    info!("omega = {omega:?}");

    // (b)
    let f = Polynomial::new([1, 1, 0, 1u128].map(R::from).to_vec());
    info!("f = {f:?}");

    let fft_f = fft(n.ilog2() as u128, f, omega);
    info!("fft(f) = {fft_f:?}");
}

#[tracing::instrument]
fn q2() {
    info!("Question 2");

    #[tracing::instrument]
    fn t(n: u64) -> u64 {
        let res = match n {
            0 => 0,
            1 => 1,
            n => t(n / 3) + t(n / 4) + n * n,
        };
        debug!("{res}");
        res
    }

    // (a)
    t(6);
}

#[tracing::instrument]
fn q3() {
    info!("Question 3");

    #[derive(Debug, Clone)]
    struct CustomOrdering;

    impl MonomialOrder<Finite<5>> for CustomOrdering {
        fn ord(
            &self,
            l: &Monomial<Finite<5>, Self>,
            r: &Monomial<Finite<5>, Self>,
        ) -> std::cmp::Ordering {
            let a = l
                .powers()
                .iter()
                .enumerate()
                .map(|(pow, n)| (pow as Natural + 1) * *n)
                .fold(0, |a, b| a + b);
            let b = r
                .powers()
                .iter()
                .enumerate()
                .map(|(pow, n)| (pow as Natural + 1) * *n)
                .fold(0, |a, b| a + b);

            b.cmp(&a).then_with(|| plex(None, l, r))
        }
    }

    // (a)

    // (b)

    let [x, y, z] = MultivariatePolynomial::init(CustomOrdering);
    let g1 = x(3) - z(1);
    let g2 = -x(1) + y(1);

    debug!("g1 = {g1:?}");
    debug!("g2 = {g2:?}");

    info!("lt(g1) = {:?}", g1.leading_term());
    info!("lt(g2) = {:?}", g2.leading_term());

    // (c)

    debug!("S(g1, g2) = {:?}", g1.s_polynomial(&g2));

    let res = multivariate_division_with_remainder(&g1.s_polynomial(&g2), &[g1, g2]);
    println!("{}", res.ascii_table());

    // (d)

    let g1 = x(4) - y(2);
    let g2 = x(1) + y(1) + z(1);

    debug!("S(g1, g2) = {:?}", g1.s_polynomial(&g2));

    let res =
        multivariate_division_with_remainder(&g1.s_polynomial(&g2), &[g1.clone(), g2.clone()]);
    println!("{}", res.ascii_table());

    let res =
        multivariate_division_with_remainder(&g1.s_polynomial(&y(2)), &[g1.clone(), g2.clone()]);
    println!("S(g1, Y^2) =>\n{}", res.ascii_table());
    let res =
        multivariate_division_with_remainder(&g2.s_polynomial(&y(2)), &[g1.clone(), g2.clone()]);
    println!("S(g2, Y^2) =>\n{}", res.ascii_table());
    let res = multivariate_division_with_remainder(&y(2), &[g1.clone(), g2.clone()]);
    println!("Y^2 =>\n{}", res.ascii_table());

    let mut b = buchbergers_algorithm(&[y(2), g1, g2]);
    println!("{b:?}");
    minimize_groebner_basis(&mut b);
    println!("{b:?}");
}

fn question_4() {}
