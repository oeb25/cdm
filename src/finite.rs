use crate::{
    field::Field,
    group::AbelianGroup,
    identity::{Addition, Identity, Multiplication},
    Group, Ring,
};

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Finite<const N: u128> {
    val: u128,
}
impl<const N: u128> std::fmt::Debug for Finite<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.val.fmt(f)
    }
}
impl<const N: u128> From<u128> for Finite<N> {
    fn from(i: u128) -> Self {
        Finite { val: i % N }
    }
}
impl<const N: u128> From<i128> for Finite<N> {
    fn from(i: i128) -> Self {
        if i >= 0 {
            Finite {
                val: (i as u128) % N,
            }
        } else {
            -Finite {
                val: i.unsigned_abs() % N,
            }
        }
    }
}

impl<const N: u128> std::ops::Add for Finite<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        ((self.val + rhs.val) % N).into()
    }
}
impl<const N: u128> std::ops::Neg for Finite<N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        ((N - self.val) % N).into()
    }
}
impl<const N: u128> std::ops::Sub for Finite<N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}
impl<const N: u128> Identity<Addition> for Finite<N> {
    fn identity() -> Self {
        Finite { val: 0 }
    }
}
impl<const N: u128> Group for Finite<N> {}
impl<const N: u128> AbelianGroup for Finite<N> {}

impl<const N: u128> std::ops::Mul for Finite<N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        ((self.val * rhs.val) % N).into()
    }
}
impl<const N: u128> Identity<Multiplication> for Finite<N> {
    fn identity() -> Self {
        Finite { val: 1 }
    }
}
impl<const N: u128> Ring for Finite<N> {
    fn multiplicative_inverse(&self) -> Option<Self> {
        (1..N)
            .find(|x| (x * self.val) % N == 1)
            .map(|x| Finite { val: x })
    }

    fn pow(&self, pow: crate::Natural) -> Self
    where
        Self: Clone,
    {
        self.val.pow(u128::from(pow) as u32).into()
    }
}

impl<const N: u128> std::ops::Div for Finite<N> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        for i in 0..N {
            let i = Finite::from(i);
            if self == rhs * i {
                return i;
            }
        }

        0u128.into()
    }
}
impl<const N: u128> Field for Finite<N> {}
