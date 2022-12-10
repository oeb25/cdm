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
        CountOps { value }
    }
}

thread_local! {
    static ADDITIONS: std::cell::Cell<u32> = Default::default();
    static MULTIPLICATIONS: std::cell::Cell<u32> = Default::default();
}

pub fn reset() {
    ADDITIONS.with(|it| it.set(0));
    MULTIPLICATIONS.with(|it| it.set(0));
}
#[derive(Debug)]
pub struct Counts {
    pub additions: u32,
    pub multiplications: u32,
}
pub fn get_counts() -> Counts {
    Counts {
        additions: ADDITIONS.with(|it| it.get()),
        multiplications: MULTIPLICATIONS.with(|it| it.get()),
    }
}

impl<F> std::ops::Add for CountOps<F>
where
    F: std::ops::Add<Output = F>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        // println!("addition");
        ADDITIONS.with(|it| it.set(it.get() + 1));
        CountOps {
            value: self.value + rhs.value,
        }
    }
}
impl<F> std::ops::Sub for CountOps<F>
where
    F: std::ops::Sub<Output = F>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        // println!("subtraction");
        ADDITIONS.with(|it| it.set(it.get() + 1));
        CountOps {
            value: self.value - rhs.value,
        }
    }
}
impl<F> std::ops::Mul for CountOps<F>
where
    F: std::ops::Mul<Output = F>,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        // println!("multiplication");
        MULTIPLICATIONS.with(|it| it.set(it.get() + 1));
        CountOps {
            value: self.value * rhs.value,
        }
    }
}
impl<F> std::ops::Div for CountOps<F>
where
    F: std::ops::Div<Output = F>,
{
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        // println!("division");
        MULTIPLICATIONS.with(|it| it.set(it.get() + 1));
        CountOps {
            value: self.value / rhs.value,
        }
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
        CountOps {
            value: F::identity(),
        }
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
        CountOps {
            value: F::identity(),
        }
    }
}
impl<F> Ring for CountOps<F>
where
    Self: AbelianGroup + std::ops::Mul<Output = Self> + Identity<Multiplication>,
    F: Ring,
{
    fn multiplicative_inverse(&self) -> Option<Self> {
        Some(CountOps {
            value: self.value.multiplicative_inverse()?,
        })
    }
}
impl<F> Field for CountOps<F> where Self: Ring + std::ops::Div<Output = Self> {}
