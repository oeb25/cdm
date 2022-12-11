use tracing::info;

fn main() {
    cdm::init_tracing();

    demo_cra();
}

fn demo_cra() {
    let res = cdm::ch05::chinese_remainder_algorithm(&[5, 7], &[1, 3]);
    info!("res = {res}");
    println!();
    println!();
    println!();
    println!();
    let res = cdm::ch05::chinese_remainder_algorithm(&[3, 4, 5], &[2, 3, 1]);
    info!("res = {res}");
}
