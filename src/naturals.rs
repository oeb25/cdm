use crate::identity::{Addition, Identity, Multiplication};

pub type Natural = u128;

impl Identity<Addition> for u128 {
    fn identity() -> Self {
        0
    }
}
impl Identity<Multiplication> for u128 {
    fn identity() -> Self {
        1
    }
}
