use cdm::{
    ch09::{fast_division_with_remainder, inversion_newton_iteration},
    Finite, Polynomial,
};
use tracing::info;

fn main() {
    cdm::init_tracing();

    demo_inversion_newton_iteration();
    demo_fast_division_with_remainder();
}

#[allow(unused)]
fn demo_inversion_newton_iteration() {
    let res = inversion_newton_iteration(
        Polynomial::<Finite<7>>::new([1, 2, 3i128].map(Finite::from).to_vec()),
        4,
    );
    info!("{res:?}");
}

fn demo_fast_division_with_remainder() {
    let res = fast_division_with_remainder(
        Polynomial::new([-4, 3, -5, 1i128].to_vec()),
        Polynomial::new([-3, 1].to_vec()),
    );
    info!("{res:?}");
}
