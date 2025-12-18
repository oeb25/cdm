use crate::{
    domain::Domain,
    euclidean_domain::EuclideanDomain,
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

macro_rules! mk_fin {
    ($($N:ident: $n:expr),*$(,)?) => {
        $(pub type $N = $crate::finite::Finite<$n>;)*
    };
}

pub mod alias {
    mk_fin!(
        F0:   0, F1:   1, F2:   2, F3:   3, F4:   4, F5:   5, F6:   6, F7:   7, F8:   8, F9:   9,
        F10: 10, F11: 11, F12: 12, F13: 13, F14: 14, F15: 15, F16: 16, F17: 17, F18: 18, F19: 19,
        F20: 20, F21: 21, F22: 22, F23: 23, F24: 24, F25: 25, F26: 26, F27: 27, F28: 28, F29: 29,
        F30: 30, F31: 31, F32: 32, F33: 33, F34: 34, F35: 35, F36: 36, F37: 37, F38: 38, F39: 39,
        F40: 40, F41: 41, F42: 42, F43: 43, F44: 44, F45: 45, F46: 46, F47: 47, F48: 48, F49: 49,
        F50: 50, F51: 51, F52: 52, F53: 53, F54: 54, F55: 55, F56: 56, F57: 57, F58: 58, F59: 59,
        F60: 60, F61: 61, F62: 62, F63: 63, F64: 64, F65: 65, F66: 66, F67: 67, F68: 68, F69: 69,
    );
}
