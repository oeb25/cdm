use cdm::{
    ch08,
    ch21::{buchbergers_algorithm, minimize_groebner_basis, multivariate_division_with_remainder},
    dft::PrimitiveRootOfUnity,
    euclidean_domain::ExtendedEuclideanAlgorithm,
    gaussian_integers::Gaussian,
    latex::ToLatex,
    mono::PLex,
    multivariate_polynomials::MultivariatePolynomial,
    rationals::rational,
    Finite, Integer, Polynomial, Rational, Ring,
};
use itertools::Itertools;

fn main() {
    tracing_subscriber::fmt::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .without_time()
        .init();

    // test_exam();
    // let m = [5, 7];
    // let v = [1, 3];
    // let res = ch05::chinese_remainder_algorithm(&m, &v);
    // println!("{:?}", res);

    // let res = cdm::ch09::inversion_newton_iteration(
    //     Polynomial::<Finite<7>>::new([1, 2, 3i128].map(Finite::from).to_vec()),
    //     4,
    // );
    // println!("{res:?}");

    // let res = cdm::ch09::fast_division_with_remainder(
    //     Polynomial::new([-4, 3, -5, 1i128].to_vec()),
    //     Polynomial::new([-3, 1].to_vec()),
    // );
    // println!("{res:?}");

    let us = [0, 1, 2, 3].map(Rational::from);
    let vs = [1, 2, 4, 8].map(Rational::from);
    let result = cdm::ch10::fast_interpolation(&us, &vs);

    eprintln!("{result:?}");
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
        let o = omega.pow(j);

        let alpha_j = f.evaluate_at(o);
        let beta_j = g.evaluate_at(o);

        let gamma_j = alpha_j * beta_j;
        dbg!(gamma_j);
        dbg!(h.evaluate_at(o));
    }

    // (v)
    let res = ch08::fast_convolution(3, f, g, omega);
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

    let [x, y, z] = MultivariatePolynomial::<Rational, _>::init(PLex(vec![2, 0, 1]));

    fn one<R: Ring>() -> R {
        R::one()
    }

    // let ord = PLex(vec![2, 0, 1]);

    let f1 = x(2) + y(2) - one();
    let f2 = rational(2.) * x(1) * y(2) - rational(2.) * x(1) * z(1) - rational(2.) * y(1);
    let f3 = x(2) - rational(2.) * y(1) * z(1) - rational(2.) * x(1) - one();

    eprintln!("{f1:?}");
    eprintln!("{f2:?}");
    // eprintln!("{f3:?}");

    let s12 = f1.s_polynomial(&f2);

    // -2*x^2*y - 2*y^2*z + 2*x*y + 2*z
    eprintln!("s12 = {s12:?}");
    let s12 = s12 * MultivariatePolynomial::one();
    eprintln!("s12 = {s12:?}");

    let fs = [f1, f2, f3];

    let table = multivariate_division_with_remainder(&s12, &fs);

    // println!("{}", table.latex_table(&ord));

    let mut res = buchbergers_algorithm(&fs);

    eprintln!("One: {res:#?}");
    minimize_groebner_basis(&mut res);
    eprintln!("Two: {res:#?}");

    let g1 = z(1) - rational(9.) * y(5) + rational(1. / 2.) * y(3);
    let g2 = x(1) + rational(9. / 2.) * y(4) + rational(1. / 2.) * y(2) - one();
    let g3 = y(6) - rational(5. / 9.) * y(4);

    eprintln!("{g1:?}");
    eprintln!("{g2:?}");
    eprintln!("{g3:?}");

    let gs = [g1, g2, g3];
    let mut basis = buchbergers_algorithm(&gs);

    eprintln!("Pre:  {basis:#?}");
    minimize_groebner_basis(&mut basis);
    eprintln!("Post: {basis:#?}");

    for i in 0..gs.len() {
        for j in i + 1..gs.len() {
            let s = gs[i].s_polynomial(&gs[j]);
            let r = multivariate_division_with_remainder(&s, &gs);

            println!("{}", r.latex_table());
        }
    }

    eprintln!("{}", gs.iter().map(|g| g.to_latex()).format("\n"));
}
