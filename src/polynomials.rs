use itertools::Itertools;

use crate::{
    euclidean_domain::EuclideanDomain,
    field::Field,
    group::AbelianGroup,
    identity::{Addition, Identity, Multiplication},
    Group, Natural, Ring,
};

#[derive(Clone)]
pub struct Polynomial<F> {
    /// The coefficients of the polynomial with the first entry being the
    /// constant term, that is:
    /// ```ignore
    /// ax^2 + bx + c == [c, b, a]
    /// ```
    coefficients: Vec<F>,
}

impl<F> PartialEq for Polynomial<F>
where
    F: PartialEq + Identity<Addition>,
{
    fn eq(&self, other: &Self) -> bool {
        self.deg() == other.deg()
            && self
                .coefficients
                .iter()
                .zip(&other.coefficients)
                .all(|(a, b)| a == b)
    }
}
impl<F> Eq for Polynomial<F> where F: PartialEq + Identity<Addition> {}

impl<F> std::fmt::Debug for Polynomial<F>
where
    F: std::fmt::Debug + Identity<Addition>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.coefficients.is_empty() {
            write!(f, "0")
        } else if Identity::<Addition>::is_identity(&self.deg()) {
            write!(f, "{:?}", self.coefficients.first().unwrap())
        } else {
            write!(
                f,
                "{}",
                self.coefficients
                    .iter()
                    .enumerate()
                    .filter(|(_, c)| !c.is_identity())
                    .map(|(i, c)| match i {
                        0 => format!("{c:?}"),
                        1 => format!("{c:?}x"),
                        _ => format!("{c:?}x^{i}"),
                    })
                    .rev()
                    .format(" + ")
            )
        }
    }
}
impl<F> Polynomial<F>
where
    F: Identity<Addition>,
{
    pub fn new(coefficients: Vec<F>) -> Self {
        let mut p = Polynomial { coefficients };
        p.normalize();
        p
    }
    pub fn normalize(&mut self) {
        self.coefficients
            .truncate(u128::from(self.deg()) as usize + 1);
    }
    pub fn normalized(&self) -> Self
    where
        F: Clone,
    {
        let mut normal = self.clone();
        normal.normalize();
        normal
    }
    pub fn deg(&self) -> Natural {
        let deg = self
            .coefficients
            .iter()
            .enumerate()
            .rfind(|(_, c)| !F::is_identity(c))
            .map(|(n, _)| n)
            .unwrap_or(0);
        Natural::from(deg as u128)
    }
    pub fn is_monic(&self) -> bool
    where
        F: Identity<Multiplication> + Clone,
    {
        Identity::<Multiplication>::is_identity(&self.coef_at(self.deg()))
    }

    pub fn scale(&self, s: &F) -> Self
    where
        F: Ring,
    {
        Self {
            coefficients: self
                .coefficients
                .iter()
                .map(|c| c.clone() * s.clone())
                .collect(),
        }
    }

    pub fn rem_pow(&self, deg: Natural) -> Self
    where
        F: Clone + std::fmt::Debug,
    {
        let mut cs = self.coefficients.clone();
        cs.truncate(deg as usize);
        Self::new(cs)
    }

    pub fn div_rem(&self, rhs: &Self) -> Option<(Self, Self)>
    where
        Self: std::fmt::Debug,
        F: Ring + std::fmt::Debug,
    {
        let a = self;
        let b = rhs;

        let n = a.deg();
        let m = b.deg();

        assert!(n >= m, "n = {n:?}, m = {m:?} ({a:?}, {b:?})");

        let mut r = a.clone();
        let u = b.lc().multiplicative_inverse()?;

        let mut q = vec![];

        for i in (0..=(n - m).into()).rev() {
            if r.deg() == m + Natural::from(i) {
                q.push(r.lc() * u.clone());
                r = r - b.times_x(i).scale(q.last().unwrap());
            } else {
                q.push(F::zero());
            }
        }

        q.reverse();

        Some((Self { coefficients: q }, r))
    }

    pub fn lc(&self) -> F
    where
        F: Group,
    {
        self.coefficients
            .last()
            .cloned()
            .unwrap_or_else(|| F::zero())
    }

    pub fn times_x(&self, pow: u128) -> Self
    where
        F: Group,
    {
        let mut new = Self {
            coefficients: self.coefficients.clone(),
        };

        for _ in 0..pow {
            new.coefficients.insert(0, F::zero());
        }

        new.normalized()
    }

    pub fn evaluate_at(&self, x: impl Into<F>) -> F
    where
        F: Ring,
    {
        let x = x.into();

        self.coefficients
            .iter()
            .rev()
            .fold(F::zero(), move |sum, c| sum * x.clone() + c.clone())
    }

    pub fn x() -> Self
    where
        F: Identity<Addition> + Identity<Multiplication>,
    {
        Self::new(vec![
            <F as Identity<Addition>>::identity(),
            <F as Identity<Multiplication>>::identity(),
        ])
    }

    pub fn diff(&self) -> Self
    where
        F: Ring,
    {
        Self::new(
            (1..=self.deg().into())
                .map(|i| {
                    self.coefficients[i as usize].clone()
                        * (0..i).map(|_| F::one()).reduce(|a, b| a + b).unwrap()
                })
                .collect(),
        )
    }

    pub fn split_poly(&self, at: usize) -> (Self, Self)
    where
        F: Clone,
    {
        (
            Self {
                coefficients: self.coefficients[0..at].to_vec(),
            },
            Self {
                coefficients: self.coefficients[at..].to_vec(),
            },
        )
    }

    pub fn rev(&self, deg: Natural) -> Self
    where
        F: Ring,
    {
        assert!(self.deg() <= deg);
        let mut coefs = self
            .iter()
            .map(|(c, _)| c.clone())
            .chain(std::iter::repeat(F::zero()))
            .take((deg + 1) as _)
            .collect_vec();
        coefs.reverse();

        Self::new(coefs)
    }

    pub fn iter(&self) -> impl Iterator<Item = (&F, Natural)> {
        self.coefficients
            .iter()
            .enumerate()
            .map(|(pow, c)| (c, (pow as u128).into()))
    }

    pub fn coef_at(&self, pow: Natural) -> F
    where
        F: Clone,
    {
        self.coefficients
            .get(u128::from(pow) as usize)
            .cloned()
            .unwrap_or_else(|| <F as Identity<Addition>>::identity())
    }
}

#[cfg(test)]
mod test {
    use crate::{ch03::ExtendedEuclideanAlgorithm, Group, Polynomial, Rational};
    use proptest::prelude::*;

    prop_compose! {
        fn polynomial()(deg in 0..5usize)(cs in prop::collection::vec(0..10i128, deg))
            -> Polynomial<Rational>
        {
            Polynomial::new(cs.into_iter().map(Rational::from).collect())
        }
    }

    proptest! {
        #[test]
        fn polynomial_division(
            a in polynomial(),
            b in polynomial(),
        ) {
            prop_assume!(a.deg() >= b.deg());

            if let Some((q, r)) = a.div_rem(&b) {
                prop_assert_eq!(a, q * b + r);
            }
        }
    }

    proptest! {
        #[test]
        fn polynomial_eed(
            a in polynomial(),
            b in polynomial(),
        ) {
            prop_assume!(a.deg() >= b.deg());
            prop_assume!(!a.is_zero() && !b.is_zero());

            prop_assert!(ExtendedEuclideanAlgorithm::property(a, b));
        }
    }

    #[test]
    fn basic_div() {
        let a = Polynomial::<Rational>::new([0, 1].map(Into::into).to_vec());
        let b = Polynomial::<Rational>::new([1].map(Into::into).to_vec());

        dbg!(a.div_rem(&b));
    }
}

impl<F> std::ops::Mul for Polynomial<F>
where
    F: Ring,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        self.coefficients
            .iter()
            .enumerate()
            .map(|(i, c)| rhs.scale(c).times_x(i as _))
            .reduce(|a, b| a + b)
            .unwrap()
    }
}
impl<F> std::ops::Mul<&Self> for Polynomial<F>
where
    F: Ring,
{
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self {
        self.coefficients
            .iter()
            .enumerate()
            .map(|(i, c)| rhs.scale(c).times_x(i as _))
            .reduce(|a, b| a + b)
            .unwrap()
    }
}
impl<F> std::ops::Mul<Self> for &Polynomial<F>
where
    F: Ring,
{
    type Output = Polynomial<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        self.coefficients
            .iter()
            .enumerate()
            .map(|(i, c)| rhs.scale(c).times_x(i as _))
            .reduce(|a, b| a + b)
            .unwrap()
    }
}
impl<F> std::ops::Add for Polynomial<F>
where
    F: Group,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self {
            coefficients: self
                .coefficients
                .iter()
                .zip_longest(&rhs.coefficients)
                .map(|cs| match cs {
                    itertools::EitherOrBoth::Both(l, r) => l.clone() + r.clone(),
                    itertools::EitherOrBoth::Left(c) | itertools::EitherOrBoth::Right(c) => {
                        c.clone()
                    }
                })
                .collect(),
        }
        .normalized()
    }
}
impl<F> std::ops::Neg for Polynomial<F>
where
    F: Clone + std::ops::Neg<Output = F> + Identity<Addition>,
{
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            coefficients: self.coefficients.into_iter().map(|c| -c).collect(),
        }
        .normalized()
    }
}
impl<F> std::ops::Sub for Polynomial<F>
where
    F: Group + std::ops::Neg<Output = F>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}
impl<F> std::ops::Rem for Polynomial<F>
where
    F: Ring + std::fmt::Debug,
{
    type Output = Self;

    fn rem(self, rhs: Self) -> Self {
        if self.deg() < rhs.deg() {
            self
        } else {
            self.div_rem(&rhs).unwrap().1.normalized()
        }
    }
}
impl<F> std::ops::Div for Polynomial<F>
where
    F: Ring + std::fmt::Debug,
{
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        self.div_rem(&rhs).unwrap().0.normalized()
    }
}
impl<F> std::ops::Sub<F> for Polynomial<F>
where
    F: Ring + std::fmt::Debug,
{
    type Output = Self;

    fn sub(self, rhs: F) -> Self {
        self - Self::new(vec![rhs])
    }
}
impl<F> std::ops::Sub<&'_ F> for Polynomial<F>
where
    F: Ring + std::fmt::Debug,
{
    type Output = Self;

    fn sub(self, rhs: &'_ F) -> Self {
        self - Self::new(vec![rhs.clone()])
    }
}

impl<F: Group> Identity<Addition> for Polynomial<F> {
    fn identity() -> Self {
        Self {
            coefficients: vec![],
        }
    }
}
impl<F: Group> Group for Polynomial<F> {}
impl<F> AbelianGroup for Polynomial<F> where F: Group {}

impl<F> Identity<Multiplication> for Polynomial<F>
where
    F: Identity<Multiplication> + Identity<Addition>,
{
    fn identity() -> Self {
        Self {
            coefficients: vec![<F as Identity<Multiplication>>::identity()],
        }
    }
}
impl<F> Ring for Polynomial<F>
where
    F: Ring,
{
    fn multiplicative_inverse(&self) -> Option<Self> {
        None
    }
}

impl<F: Field + std::fmt::Debug> EuclideanDomain for Polynomial<F> {
    fn d(&self) -> Option<Natural> {
        Some(self.deg())
    }
}
