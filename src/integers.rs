use crate::prelude::*;

pub type Integer = i128;

impl Identity<Addition> for i128 {
    fn identity() -> Self {
        0
    }
}
impl Group for i128 {}
impl AbelianGroup for i128 {}
impl Identity<Multiplication> for i128 {
    fn identity() -> Self {
        1
    }
}
impl Ring for i128 {
    fn multiplicative_inverse(&self) -> Option<Self> {
        if *self == -1 {
            Some(-1)
        } else if *self == 1 {
            Some(1)
        } else {
            None
        }
    }
}
impl Domain for i128 {}
impl EuclideanDomain for i128 {
    fn d(&self) -> Option<Natural> {
        Some(self.unsigned_abs())
    }
}
