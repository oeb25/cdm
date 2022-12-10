use crate::ring::Ring;

/// Let (R, +, ·) be a commutative rings such that R∗ = R\{0}. Then R is called a field.
pub trait Field: Ring + std::ops::Div<Output = Self> {}
