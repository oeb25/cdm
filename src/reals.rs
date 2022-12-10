use derive_more::{Add, Div, Mul, Neg, Sub};

use crate::{
    field::Field,
    group::AbelianGroup,
    identity::{Addition, Identity, Multiplication},
    Group, Ring,
};

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Add, Sub, Neg, Mul, Div)]
#[mul(forward)]
#[div(forward)]
pub struct Real {
    value: f64,
}
impl Real {
    pub fn abs(self) -> Self {
        Real {
            value: self.value.abs(),
        }
    }
}

impl Identity<Addition> for Real {
    fn identity() -> Self {
        Real { value: 0.0 }
    }
}
impl Group for Real {}
impl AbelianGroup for Real {}
impl Identity<Multiplication> for Real {
    fn identity() -> Self {
        Real { value: 1.0 }
    }
}
impl Ring for Real {
    fn multiplicative_inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(Real {
                value: 1.0 / self.value,
            })
        }
    }
}
impl Field for Real {}
