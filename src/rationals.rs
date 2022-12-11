use crate::{
    field::Field,
    group::AbelianGroup,
    identity::{Addition, Identity, Multiplication},
    Group, Integer, Natural, Ring,
};

pub fn rational(f: f64) -> Rational {
    Rational::approximate(f)
}

#[derive(Clone, Copy, Eq, PartialOrd)]
pub struct Rational {
    pub num: Integer,
    pub denom: Natural,
}

impl PartialEq for Rational {
    fn eq(&self, other: &Self) -> bool {
        i128::from(self.num) as f64 / u128::from(self.denom) as f64
            == i128::from(other.num) as f64 / u128::from(other.denom) as f64
    }
}
impl std::hash::Hash for Rational {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        (((i128::from(self.num) as f64 / u128::from(self.denom) as f64) * 1000000.0) as i128)
            .hash(state);
    }
}

impl std::fmt::Debug for Rational {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if u128::from(self.denom) == 1 {
            write!(f, "{:?}", self.num)
        } else {
            write!(f, "{:?}/{:?}", self.num, self.denom)
        }
    }
}
impl Rational {
    pub fn approximate(f: f64) -> Self {
        let g = f.abs();

        let mut num = 0;
        let mut denom = 1;

        while (num as f64 / denom as f64 - g).abs() > f64::EPSILON {
            let delta = num as f64 / denom as f64 - g;

            if delta < 0.0 {
                num += 1;
            } else {
                denom += 1;
            }
        }

        if f < 0.0 {
            num = -num;
        }

        Rational {
            num: num as _,
            denom,
        }
    }
    pub fn abs(self) -> Self {
        Rational {
            num: self.num.abs(),
            denom: self.denom,
        }
    }
    pub fn normalized(self) -> Self {
        if self.num.is_zero() {
            return Self::zero();
        }

        let mut n = u128::from(self.num.unsigned_abs());
        let mut m = u128::from(self.denom);

        while m != 0 {
            let tmp = n;
            n = m;
            m = tmp % m;
        }

        let gcd = n;

        Rational {
            num: self.num / Integer::from(gcd as i128),
            denom: self.denom / Natural::from(gcd),
        }
    }
}

#[test]
fn rational_signed_normalize() {
    let minus_1 = Rational {
        num: (-1i128).into(),
        denom: 1u128.into(),
    }
    .normalized();

    todo!("{:?}", minus_1 / minus_1);
}

impl From<i128> for Rational {
    fn from(n: i128) -> Self {
        Rational {
            num: n.into(),
            denom: 1u128.into(),
        }
    }
}

impl std::ops::Add for Rational {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Rational {
            num: self.num * rhs.denom as i128 + self.denom as i128 * rhs.num,
            denom: self.denom * rhs.denom,
        }
        .normalized()
    }
}
impl std::ops::Sub for Rational {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Rational {
            num: self.num * rhs.denom as i128 - self.denom as i128 * rhs.num,
            denom: self.denom * rhs.denom,
        }
        .normalized()
    }
}
impl std::ops::Neg for Rational {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Rational {
            num: -self.num,
            denom: self.denom,
        }
        .normalized()
    }
}
impl std::ops::Rem for Rational {
    type Output = Self;

    fn rem(self, _: Self) -> Self::Output {
        Self::zero()
    }
}
impl std::ops::Div for Rational {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Rational {
            num: if rhs.num == rhs.num.abs() {
                self.num * rhs.denom as i128
            } else {
                -self.num * rhs.denom as i128
            },
            denom: self.denom * rhs.num.unsigned_abs(),
        }
        .normalized()
    }
}
impl std::ops::Mul for Rational {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Rational {
            num: self.num * rhs.num,
            denom: self.denom * rhs.denom,
        }
        .normalized()
    }
}

impl Identity<Addition> for Rational {
    fn identity() -> Self {
        Rational {
            num: 0i128.into(),
            denom: 1u128.into(),
        }
    }
}
impl Group for Rational {}
impl AbelianGroup for Rational {}

impl Identity<Multiplication> for Rational {
    fn identity() -> Self {
        Rational {
            num: 1i128.into(),
            denom: 1u128.into(),
        }
    }
}
impl Ring for Rational {
    fn multiplicative_inverse(&self) -> Option<Self> {
        Some(match i128::from(self.num) {
            0 => 0.into(),
            x if x < 0 => Rational {
                num: -(self.denom as i128),
                denom: self.num.unsigned_abs(),
            },
            _ => Rational {
                num: self.denom as i128,
                denom: self.num.unsigned_abs(),
            },
        })
    }
}
impl Field for Rational {}
