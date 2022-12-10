use derive_more::{Add, Div, From, Into, Mul, Neg, Rem, Sub};

use crate::{
    domain::Domain,
    euclidean_domain::EuclideanDomain,
    group::AbelianGroup,
    identity::{Addition, Identity, Multiplication},
    Group, Integer, Natural, Ring,
};

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Hash, Add, Sub, Neg, From, Into)]
// #[mul(forward)]
// #[div(forward)]
// #[rem(forward)]
pub struct Gaussian<F> {
    pub a: F,
    pub b: F,
}
impl<F: std::fmt::Debug> std::fmt::Debug for Gaussian<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?} + {:?}i", self.a, self.b)
    }
}

impl<F> Gaussian<F> {
    pub fn new(a: F, b: F) -> Self {
        Gaussian { a, b }
    }
    pub fn map<T>(self, f: impl Fn(F) -> T) -> Gaussian<T> {
        Gaussian {
            a: f(self.a),
            b: f(self.b),
        }
    }
    pub fn conj(self) -> Self
    where
        F: std::ops::Neg<Output = F>,
    {
        Gaussian {
            a: self.a,
            b: -self.b,
        }
    }
}

impl<F> std::ops::Mul for Gaussian<F>
where
    F: std::ops::Add<Output = F> + std::ops::Sub<Output = F> + std::ops::Mul<Output = F> + Clone,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        // (a + bi)(c + di)
        // (ac + adi + bic + bidi)
        // (ac - bd + (ad + bc)i)

        Gaussian {
            a: self.a.clone() * rhs.a.clone() - self.b.clone() * rhs.b.clone(),
            b: self.a.clone() * rhs.b + self.b * rhs.a,
        }
    }
}
impl<F> std::ops::Div for Gaussian<F>
where
    F: std::fmt::Debug
        + Clone
        + std::ops::Add<Output = F>
        + std::ops::Sub<Output = F>
        + std::ops::Mul<Output = F>
        + std::ops::Div<Output = F>,
    Self: EuclideanDomain,
{
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        // (a + bi)/(c + di) = (ac + bd + (bc - ad)i)/(cc + dd)
        // ((a + bi)*(c - di))/((c + di)*(c - di))
        // ((a + bi)*(c - di))/(cc + dd)
        // (ac - adi + bic - bidi)/(cc + dd)
        // (ac - adi + bic + bd)/(cc + dd)
        // (ac + bd + (bc - ad)i)/(cc + dd)

        let a = self.a.clone();
        let b = self.b.clone();
        let c = rhs.a.clone();
        let d = rhs.b.clone();

        let denum = c.clone() * c.clone() + d.clone() * d.clone();

        let res = Gaussian {
            a: (a.clone() * c.clone() + b.clone() * d.clone()) / denum.clone(),
            b: (b * c - a * d) / denum,
        };

        println!("{self:?} / {rhs:?}");
        assert!(res.d() < self.d(), "{res:?} !< {self:?}");

        res
    }
}
impl<F> std::ops::Rem for Gaussian<F> {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        todo!()
    }
}

impl<F> Identity<Addition> for Gaussian<F>
where
    F: Identity<Addition>,
{
    fn identity() -> Self {
        Gaussian {
            a: F::identity(),
            b: F::identity(),
        }
    }
}
impl<F: Group> Group for Gaussian<F> {}
impl<F: Group> AbelianGroup for Gaussian<F> {}
impl<F: Identity<Multiplication> + Identity<Addition>> Identity<Multiplication> for Gaussian<F> {
    fn identity() -> Self {
        Gaussian {
            a: <F as Identity<Multiplication>>::identity(),
            b: <F as Identity<Addition>>::identity(),
        }
    }
}
impl<F: Ring> Ring for Gaussian<F> {
    fn multiplicative_inverse(&self) -> Option<Self> {
        // TODO
        // None
        todo!()
    }
}
impl<F: Domain> Domain for Gaussian<F> {}
impl EuclideanDomain for Gaussian<Integer> {
    fn d(&self) -> Option<Natural> {
        Some(self.a.abs_nat().pow(2u128.into()) + self.b.abs_nat().pow(2u128.into()))
    }
}
