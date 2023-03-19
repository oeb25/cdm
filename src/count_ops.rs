use derive_more::Neg;

use crate::{
    field::Field,
    group::AbelianGroup,
    identity::{Addition, Identity, Multiplication},
    Group, Ring,
};

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Neg, Hash)]
pub struct CountOps<F> {
    pub value: F,
}
impl<F> From<F> for CountOps<F> {
    fn from(value: F) -> Self {
        Self { value }
    }
}

thread_local! {
    static COUNTS:  std::cell::Cell<Counts> = Default::default();
}

pub fn reset() {
    COUNTS.with(|x| x.set(Default::default()))
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Counts {
    pub add: u32,
    pub mul: u32,
}

impl Counts {
    fn inc_add(&self) -> Self {
        Self {
            add: self.add + 1,
            mul: self.mul,
        }
    }

    fn inc_mul(&self) -> Self {
        Self {
            add: self.add,
            mul: self.mul + 1,
        }
    }
}

pub fn inc_add() {
    COUNTS.with(|x| x.get().inc_add());
}

pub fn inc_mul() {
    COUNTS.with(|x| x.get().inc_mul());
}

pub fn get_counts() -> Counts {
    COUNTS.with(|x| x.get())
}

impl<F> std::ops::Add for CountOps<F>
where
    F: std::ops::Add<Output = F>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        inc_add();
        (self.value + rhs.value).into()
    }
}
impl<F> std::ops::Sub for CountOps<F>
where
    F: std::ops::Sub<Output = F>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        inc_add();
        (self.value - rhs.value).into()
    }
}
impl<F> std::ops::Mul for CountOps<F>
where
    F: std::ops::Mul<Output = F>,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        inc_mul();
        (self.value * rhs.value).into()
    }
}
impl<F> std::ops::Div for CountOps<F>
where
    F: std::ops::Div<Output = F>,
{
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        inc_mul();
        (self.value / rhs.value).into()
    }
}

impl<F> std::fmt::Debug for CountOps<F>
where
    F: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.value.fmt(f)
    }
}

// impl<F> std::ops::Add for CountOps<F> {
//     type Output = Self;

//     fn add(self, rhs: Self) -> Self::Output {
//         CountOps {
//             value: self.value + rhs.value,
//         }
//     }
// }
// impl<F> std::ops::Sub for CountOps<F> {
//     type Output = Self;

//     fn sub(self, rhs: Self) -> Self::Output {
//         CountOps {
//             value: self.value - rhs.value,
//         }
//     }
// }
// impl<F> std::ops::Neg for CountOps<F> {
//     type Output = Self;

//     fn neg(self) -> Self::Output {
//         CountOps { value: -self.value }
//     }
// }
// impl<F> std::ops::Rem for CountOps<F> {
//     type Output = Self;

//     fn rem(self, _: Self) -> Self::Output {
//         Self::zero()
//     }
// }
// impl<F> std::ops::Div for CountOps<F> {
//     type Output = Self;

//     fn div(self, rhs: Self) -> Self::Output {
//         CountOps {
//             value: self.value - rhs.value,
//         }
//     }
// }
// impl<F> std::ops::Mul for CountOps<F> {
//     type Output = Self;

//     fn mul(self, rhs: Self) -> Self::Output {
//         CountOps {
//             num: self.num * rhs.num,
//             denom: self.denom * rhs.denom,
//         }
//         .normalized()
//     }
// }

impl<F> Identity<Addition> for CountOps<F>
where
    F: Identity<Addition>,
{
    fn identity() -> Self {
        F::identity().into()
    }
}
impl<F> Group for CountOps<F> where
    Self: std::ops::Add<Self, Output = Self>
        + std::ops::Neg<Output = Self>
        + std::ops::Sub<Output = Self>
        + Identity<Addition>
        + Sized
        + Clone
        + PartialEq
        + std::fmt::Debug
{
}
impl<F> AbelianGroup for CountOps<F>
where
    Self: Group,
    F: AbelianGroup,
{
}

impl<F> Identity<Multiplication> for CountOps<F>
where
    F: Identity<Multiplication>,
{
    fn identity() -> Self {
        F::identity().into()
    }
}
impl<F> Ring for CountOps<F>
where
    Self: AbelianGroup + std::ops::Mul<Output = Self> + Identity<Multiplication>,
    F: Ring,
{
    fn multiplicative_inverse(&self) -> Option<Self> {
        Self {
            value: self.value.multiplicative_inverse()?,
        }
        .into()
    }
}
impl<F> Field for CountOps<F> where Self: Ring + std::ops::Div<Output = Self> {}
