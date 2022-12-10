use crate::identity::{Addition, Identity};

pub trait Group:
    std::ops::Add<Self, Output = Self>
    + std::ops::Neg<Output = Self>
    + std::ops::Sub<Output = Self>
    + Identity<Addition>
    + Sized
    + Clone
    + PartialEq
    + std::fmt::Debug
{
    /// The identity element for addition
    fn zero() -> Self {
        Self::identity()
    }
    fn is_zero(&self) -> bool {
        &Self::zero() == self
    }
}

/// Is a group for which the add operation is commutative
pub trait AbelianGroup: Group {}
