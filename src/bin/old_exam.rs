use cdm::{
    dft::{fft, PrimitiveRootOfUnity},
    groebner,
    mono::{MonomialOrder, PLex},
    mpoly,
    multivariate_polynomials::{self, MultivariatePolynomial},
    Finite, Monomial, Natural, Polynomial,
};
use tracing::{debug, info};

fn main() {
    tracing_subscriber::fmt::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .without_time()
        .init();

    question_1();
    question_2();
    question_3();
}

fn question_1() {
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

fn question_2() {
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

fn question_3() {
    struct CustomOrdering;

    impl MonomialOrder<Finite<5>> for CustomOrdering {
        fn ord(&self, l: &Monomial<Finite<5>>, r: &Monomial<Finite<5>>) -> std::cmp::Ordering {
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

            b.cmp(&a).then_with(|| PLex::default().ord(l, r))
        }
    }

    // (a)

    // (b)

    type F = Finite<5>;

    let g1 = mpoly::<F, 2>([(1i128, &[3]), (-1i128, &[0, 0, 1])]);
    let g2 = mpoly::<F, 2>([(-1i128, &[1]), (1i128, &[0, 1])]);

    debug!("g1 = {g1:?}");
    debug!("g2 = {g2:?}");

    info!("lt(g1) = {:?}", g1.leading_term(&CustomOrdering));
    info!("lt(g2) = {:?}", g2.leading_term(&CustomOrdering));

    // (c)

    let res = multivariate_polynomials::multivariate_division_with_remainder(
        &CustomOrdering,
        &g1.s_polynomial(&g2, &CustomOrdering),
        &[g1, g2],
    );

    println!("{}", res.ascii_table(&CustomOrdering));
}

fn question_4() {}
