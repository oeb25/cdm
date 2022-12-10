use derive_more::{Add, Div, From, Into, Mul, Neg, Rem, Sub};

use crate::{
    euclidean_domain::EuclideanDomain,
    group::AbelianGroup,
    identity::{Addition, Identity, Multiplication},
    Group, Natural, Ring,
};

#[derive(
    Clone, Copy, PartialEq, Eq, PartialOrd, Hash, Add, Sub, Neg, Rem, Mul, Div, From, Into,
)]
#[mul(forward)]
#[div(forward)]
#[rem(forward)]
pub struct Integer {
    value: i128,
}
impl std::fmt::Debug for Integer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.value.fmt(f)
    }
}
impl Integer {
    pub fn abs(self) -> Self {
        Integer {
            value: self.value.abs(),
        }
    }
    pub fn abs_nat(self) -> Natural {
        self.value.unsigned_abs().into()
    }
}
impl From<Natural> for Integer {
    fn from(n: Natural) -> Self {
        Integer {
            value: u128::from(n) as _,
        }
    }
}

impl Identity<Addition> for Integer {
    fn identity() -> Self {
        Integer { value: 0 }
    }
}
impl Group for Integer {}
impl AbelianGroup for Integer {}
impl Identity<Multiplication> for Integer {
    fn identity() -> Self {
        1.into()
    }
}
impl Ring for Integer {
    fn multiplicative_inverse(&self) -> Option<Self> {
        if self.value == -1 {
            Some(Integer { value: -1 })
        } else if self.value == 1 {
            Some(Integer { value: 1 })
        } else {
            None
        }
    }
}
impl EuclideanDomain for Integer {
    fn d(&self) -> Option<Natural> {
        Some(self.abs_nat())
    }
}
