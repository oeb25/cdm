use crate::{
    field::Field,
    group::AbelianGroup,
    identity::{Addition, Identity, Multiplication},
    Group, Integer, Natural, Ring,
};

pub fn finite<const N: Natural>(x: Integer) -> Finite<N> {
    Finite::from(x)
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Finite<const N: Natural> {
    val: Natural,
}
impl<const N: Natural> std::fmt::Debug for Finite<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.val.fmt(f)
    }
}
impl<const N: Natural> From<Natural> for Finite<N> {
    fn from(i: Natural) -> Self {
        Finite { val: i % N }
    }
}
impl<const N: Natural> From<Integer> for Finite<N> {
    fn from(i: Integer) -> Self {
        if i >= 0 {
            Finite {
                val: (i as Natural) % N,
            }
        } else {
            -Finite {
                val: i.unsigned_abs() % N,
            }
        }
    }
}

impl<const N: Natural> std::ops::Add for Finite<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        ((self.val + rhs.val) % N).into()
    }
}
impl<const N: Natural> std::ops::Neg for Finite<N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        ((N - self.val) % N).into()
    }
}
impl<const N: Natural> std::ops::Sub for Finite<N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}
impl<const N: Natural> Identity<Addition> for Finite<N> {
    fn identity() -> Self {
        Finite { val: 0 }
    }
}
impl<const N: Natural> Group for Finite<N> {}
impl<const N: Natural> AbelianGroup for Finite<N> {}

impl<const N: Natural> std::ops::Mul for Finite<N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        ((self.val * rhs.val) % N).into()
    }
}
impl<const N: Natural> Identity<Multiplication> for Finite<N> {
    fn identity() -> Self {
        Finite { val: 1 }
    }
}
impl<const N: Natural> Ring for Finite<N> {
    fn multiplicative_inverse(&self) -> Option<Self> {
        (1..N)
            .find(|x| (x * self.val) % N == 1)
            .map(|x| Finite { val: x })
    }

    fn pow(&self, pow: crate::Natural) -> Self
    where
        Self: Clone,
    {
        self.val.pow(Natural::from(pow) as u32).into()
    }
}

impl<const N: Natural> std::ops::Div for Finite<N> {
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
impl<const N: Natural> Field for Finite<N> {}
