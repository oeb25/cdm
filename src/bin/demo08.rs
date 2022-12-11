use cdm::{
    ch08::{
        fast_convolution, fast_fourier_transform, karatsubas_polynomial_multiplication_algorithm,
    },
    dft::PrimitiveRootOfUnity,
    Finite, Integer, Polynomial, Ring,
};
use tracing::info;

fn main() {
    cdm::init_tracing();

    // demo_fast_fourier_transform();
    // demo_fast_convolution();
    demo_karatsubas_polynomial_multiplication_algorithm();
}

#[allow(unused)]
// cdm::ch08::karatsubas_polynomial_multiplication_algorithm(k, f, g);
fn demo_karatsubas_polynomial_multiplication_algorithm() {
    type R = Integer;

    let f = Polynomial::new([3, -4, 3, 5i128].map(R::from).to_vec());
    let g = Polynomial::new([-2, 7, -5, 2i128].map(R::from).to_vec());

    karatsubas_polynomial_multiplication_algorithm(2, f, g);
}

#[allow(unused)]
/// Old exam 2021 Question 1
fn demo_fast_fourier_transform() {
    type R = Finite<5>;
    let n = 4;

    // (a)
    let omega = PrimitiveRootOfUnity::new(n, R::from(2i128)).expect("was not a root of unity");
    info!("omega = {omega:?}");
    info!("omega^-1 = {:?}", omega.inverse());

    // (b)
    let f = Polynomial::new([1, 1, 0, 1u128].map(R::from).to_vec());
    info!("f = {f:?}");

    let fft_f = fast_fourier_transform(n.ilog2() as u128, f, omega);
    info!("fft(f) = {fft_f:?}");
}

#[allow(unused)]
/// Exercise 8.10
fn demo_fast_convolution() {
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
        // debug!("gamma_{j} = {:?}", gamma_j);
        // debug!("h({o:?}) = {:?}", h.evaluate_at(o));
    }

    // (v)
    let res = fast_convolution(3, f, g, omega);
    info!("h   = {:?}", h);
    info!("res = {:?}", res);
}

// cdm::ch08::fast_negative_wrapped_convolution();
