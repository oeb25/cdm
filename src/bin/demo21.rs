use cdm::{
    ch21::{buchbergers_algorithm, multivariate_division_with_remainder},
    mono::PLex,
    multivariate_polynomials::MultivariatePolynomial,
    Finite,
};
use tracing::{debug, info};

fn main() {
    cdm::init_tracing();

    demo_leading_term();
    demo_s_polynomial();
    demo_multivariant_div_w_rem();
    demo_buchbergers_algorithm();
    // cdm::ch21::multivariate_division_with_remainder(f, fs);
}

fn demo_leading_term() {
    let [x, y, z] = MultivariatePolynomial::<Finite<5>, _>::init(PLex::default());
    let g1 = x(3) - z(1);
    let g2 = -x(1) + y(1);

    debug!("g1 = {g1:?}");
    debug!("g2 = {g2:?}");

    info!("lt(g1) = {:?}", g1.leading_term());
    info!("lt(g2) = {:?}", g2.leading_term());
}

fn demo_s_polynomial() {
    let [x, y, z] = MultivariatePolynomial::<Finite<5>, _>::init(PLex::default());
    let g1 = x(3) - z(1);
    let g2 = -x(1) + y(1);
    info!("S(g1, g2) = {:?}", g1.s_polynomial(&g2));
}

fn demo_multivariant_div_w_rem() {
    let [x, y, z] = MultivariatePolynomial::<Finite<5>, _>::init(PLex::default());
    let g1 = x(3) - z(1);
    let g2 = -x(1) + y(1);
    let res = multivariate_division_with_remainder(&g1.s_polynomial(&g2), &[g1, g2]);
    println!("{}", res.ascii_table());
}

fn demo_buchbergers_algorithm() {
    let [x, y, z] = MultivariatePolynomial::<Finite<5>, _>::init(PLex::default());
    let g1 = x(3) - z(1);
    let g2 = -x(1) + y(1);
    let res = buchbergers_algorithm(&[g1, g2]);
    debug!("{:?}", res);
}
