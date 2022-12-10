use derive_more::{Add, Div, From, Into, Mul, Rem, Sub};

use crate::identity::{Addition, Identity, Multiplication};

#[derive(
    Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Add, Sub, Rem, Mul, Div, From, Into,
)]
#[mul(forward)]
#[div(forward)]
#[rem(forward)]
pub struct Natural {
    value: u128,
}
impl Natural {
    pub const fn zero() -> Self {
        Natural { value: 0 }
    }
    pub const fn is_zero(self) -> bool {
        self.value == 0
    }

    pub fn pow(&self, k: Self) -> Self {
        Natural {
            value: self.value.pow(k.value as _),
        }
    }
}
impl std::fmt::Debug for Natural {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.value.fmt(f)
    }
}
impl From<Natural> for u64 {
    fn from(val: Natural) -> Self {
        val.value as _
    }
}
impl From<Natural> for usize {
    fn from(val: Natural) -> Self {
        val.value as _
    }
}
impl From<Natural> for i128 {
    fn from(val: Natural) -> Self {
        val.value as _
    }
}

impl Identity<Addition> for Natural {
    fn identity() -> Self {
        Natural { value: 0 }
    }
}
impl Identity<Multiplication> for Natural {
    fn identity() -> Self {
        Natural { value: 1 }
    }
}
