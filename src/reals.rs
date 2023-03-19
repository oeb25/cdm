use derive_more::{Add, Div, Mul, Neg, Sub};

use crate::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Add, Sub, Neg, Mul, Div)]
#[mul(forward)]
#[div(forward)]
pub struct Real(f64);
impl Real {
    pub fn abs(self) -> Self {
        Self(self.0.abs())
    }
}

impl Identity<Addition> for Real {
    fn identity() -> Self {
        Self(0.0)
    }
}
impl Group for Real {}
impl AbelianGroup for Real {}
impl Identity<Multiplication> for Real {
    fn identity() -> Self {
        Self(1.0)
    }
}
impl Ring for Real {
    fn multiplicative_inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Self(1.0 / self.0).into()
        }
    }
}
impl Field for Real {}
