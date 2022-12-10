#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Addition;
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Multiplication;

pub trait Identity<Op>: PartialEq {
    fn identity() -> Self;
    fn is_identity(&self) -> bool
    where
        Self: PartialEq + Sized,
    {
        self == &Self::identity()
    }
}
