use cdm::{
    ch10::{fast_interpolation, fast_multipoint_evaluation},
    Polynomial, Rational,
};
use itertools::Itertools;
use tracing::info;

fn main() {
    cdm::init_tracing();

    demo_fast_multipoint_evaluation();
    demo_fast_interpolation();
}

fn demo_fast_multipoint_evaluation() {
    let f = Polynomial::new([1, 2, 3].map(Rational::from).to_vec());
    let us = (-3..=4).map(Rational::from).collect_vec();

    let eval = fast_multipoint_evaluation(f, &us);

    info!("{eval:?}");
}

fn demo_fast_interpolation() {
    let us = [0, 1, 2, 3].map(Rational::from);
    let vs = [1, 2, 4, 8].map(Rational::from);
    let result = fast_interpolation(&us, &vs);

    info!("{result:?}");
}
