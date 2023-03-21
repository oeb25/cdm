use either::Either;
use itertools::Itertools;

use crate::prelude::*;

#[derive(Clone)]
pub enum MultivariatePolynomial<F, O: MonomialOrder<F>> {
    Constant(F),
    Terms { ord: O, terms: Vec<Monomial<F, O>> },
}

impl<F, O> PartialEq for MultivariatePolynomial<F, O>
where
    F: Ring,
    O: MonomialOrder<F>,
{
    fn eq(&self, other: &Self) -> bool {
        use itertools::EitherOrBoth;
        self.terms().zip_longest(other.terms()).all(|ps| {
            let EitherOrBoth::Both(l, r) = ps else { return false };
            l == r
        })
    }
}

impl<F, O> std::fmt::Debug for MultivariatePolynomial<F, O>
where
    F: std::fmt::Debug + Ring,
    O: MonomialOrder<F>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = format!("{:?}", self.terms().filter(|t| !t.is_zero()).format(" + "));
        if s.is_empty() {
            write!(f, "0")
        } else {
            write!(f, "{}", s)
        }
    }
}

impl<F, O> MultivariatePolynomial<F, O>
where
    F: Ring,
    O: MonomialOrder<F>,
{
    pub fn init<const N: usize>(ord: O) -> [impl Fn(Natural) -> Self; N] {
        std::array::from_fn(|i| {
            let ord = ord.clone();
            move |p| {
                let mut vs = vec![0; i + 1];
                vs[i] = p;
                Self::new(ord.clone(), vec![Monomial::new(ord.clone(), F::one(), vs)])
            }
        })
    }
    pub fn new(ord: O, mut terms: Vec<Monomial<F, O>>) -> Self {
        terms.sort_by(|l, r| ord.ord(r, l));

        let mut terms = terms
            .into_iter()
            .group_by(|t| t.without_coef())
            .into_iter()
            .map(|(a, b)| a * b.map(|t| t.coef().clone()).reduce(|a, b| a + b).unwrap())
            .collect_vec();
        terms.sort_by(|l, r| ord.ord(r, l));

        Self::Terms { ord, terms }
    }
    pub fn try_new(ord: Option<O>, terms: Vec<Monomial<F, O>>) -> Option<Self> {
        let ord = ord.or_else(|| terms.iter().find_map(|t| t.ord().cloned()));

        if let Some(ord) = ord {
            Some(Self::new(ord, terms))
        } else {
            let mut coef = F::zero();

            for t in &terms {
                if !t.powers().is_empty() || t.powers().iter().any(|p| *p != 0) {
                    return None;
                }
                coef = coef + t.coef().clone();
            }
            Some(Self::Constant(coef))
        }
    }
    pub fn ord(&self) -> Option<&O> {
        let Self::Terms { ord, ..} = self else { return None };
        ord.into()
    }

    pub fn leading_term(&self) -> Monomial<F, O> {
        self.terms()
            .max_by(|l, r| self.ord().unwrap().ord(l, r))
            .unwrap_or_else(|| Monomial::zero(self.ord().cloned()))
    }

    pub fn leading_coef(&self) -> F {
        self.leading_term().coef().clone()
    }

    pub fn leading_monomial(&self) -> Monomial<F, O>
    where
        F: Identity<Multiplication>,
    {
        self.leading_term().without_coef()
    }

    pub fn multi_deg(&self) -> Vec<Natural> {
        todo!()
    }

    pub fn terms(&self) -> impl Iterator<Item = Monomial<F, O>> + Clone + '_ {
        match self {
            Self::Constant(c) => {
                Either::Left([Monomial::constant(self.ord().cloned(), c.clone())].into_iter())
            }
            Self::Terms { terms, .. } => Either::Right(terms.iter().cloned()),
        }
        .filter(|t| !t.is_zero())
    }

    pub fn minimize(&self) -> Self
    where
        F: Field,
    {
        let terms = self
            .terms()
            .map(|t| t.clone() / self.leading_coef())
            .collect();

        Self::try_new(self.ord().cloned(), terms).unwrap()
    }

    pub fn s_polynomial(&self, other: &Self) -> Self
    where
        F: Field + std::fmt::Debug,
    {
        let ord = self.ord().or_else(|| other.ord()).unwrap();

        let lcm = self
            .leading_term()
            .powers()
            .iter()
            .copied()
            .zip_longest(other.leading_term().powers().iter().copied())
            .map(|ps| ps.reduce(|a, b| a.max(b)))
            .collect_vec();

        let lcm = Monomial::new(ord.clone(), F::one(), lcm.clone());

        let l = lcm
            .div(&self.leading_term())
            .unwrap_or_else(|| panic!("{lcm:?}/{:?} failed", self.leading_term()));
        let r = lcm
            .div(&other.leading_term())
            .unwrap_or_else(|| panic!("{lcm:?}/{:?} failed", other.leading_term()));

        let l = Self::from(l);
        let r = Self::from(r);

        l * self.clone() - r * other.clone()
    }

    pub fn div_mono(&self, m: &Monomial<F, O>) -> Option<Self>
    where
        F: Field + std::fmt::Debug,
    {
        let terms = self.terms().map(|t| t.div(m)).collect::<Option<_>>()?;

        Self::try_new(self.ord().or_else(|| m.ord()).cloned(), terms)
    }

    pub fn constant(ord: Option<O>, c: impl Into<F>) -> Self {
        let Some(ord) = ord else { return Self::Constant(c.into());};
        Self::new(ord.clone(), vec![Monomial::new(ord, c, vec![])])
    }
}

impl<F, O> From<Monomial<F, O>> for MultivariatePolynomial<F, O>
where
    F: Ring,
    O: MonomialOrder<F>,
{
    fn from(value: Monomial<F, O>) -> Self {
        Self::try_new(value.ord().cloned(), vec![value]).unwrap()
    }
}

impl<F, O> std::ops::Add for MultivariatePolynomial<F, O>
where
    F: Ring,
    O: MonomialOrder<F>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self::try_new(
            self.ord().or_else(|| rhs.ord()).cloned(),
            self.terms()
                .into_iter()
                .chain(rhs.terms().into_iter())
                .collect(),
        )
        .unwrap()
    }
}
impl<F, O> std::ops::Neg for MultivariatePolynomial<F, O>
where
    F: Ring,
    O: MonomialOrder<F>,
{
    type Output = Self;

    fn neg(self) -> Self {
        Self::try_new(
            self.ord().cloned(),
            self.terms()
                .into_iter()
                .map(|t| t.map_coef(|c| -c.clone()))
                .collect(),
        )
        .unwrap()
    }
}
impl<F, O> std::ops::Sub for MultivariatePolynomial<F, O>
where
    F: Ring,
    O: MonomialOrder<F>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}

impl<F: Ring, O> Identity<Addition> for MultivariatePolynomial<F, O>
where
    O: MonomialOrder<F>,
{
    fn identity() -> Self {
        Self::Constant(F::zero())
    }
}
impl<F: Ring, O> Group for MultivariatePolynomial<F, O> where O: MonomialOrder<F> {}
impl<F: Ring, O> AbelianGroup for MultivariatePolynomial<F, O> where O: MonomialOrder<F> {}

impl<F: Ring, O> Identity<Multiplication> for MultivariatePolynomial<F, O>
where
    O: MonomialOrder<F>,
{
    fn identity() -> Self {
        Self::Constant(F::one())
    }
}
impl<F: Ring, O> Ring for MultivariatePolynomial<F, O>
where
    O: MonomialOrder<F>,
{
    fn multiplicative_inverse(&self) -> Option<Self> {
        todo!()
    }
}

impl<F: Ring, O> std::ops::Mul for MultivariatePolynomial<F, O>
where
    O: MonomialOrder<F>,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let terms = self
            .terms()
            .cartesian_product(rhs.terms())
            .map(|(a, b)| a.clone() * b.clone())
            .collect();

        Self::try_new(self.ord().or_else(|| rhs.ord()).cloned(), terms).unwrap()
    }
}

impl<O, const N: Natural> std::ops::Mul<MultivariatePolynomial<Self, O>> for Finite<N>
where
    O: MonomialOrder<Self>,
{
    type Output = MultivariatePolynomial<Self, O>;

    fn mul(self, rhs: Self::Output) -> Self::Output {
        Self::Output::constant(None, self) * rhs
    }
}

// we implement this via a macro rather than by a blanket implementation
// like below since the implementation rules in Rust are kind of stupid
macro_rules! impl_mul {
    ($ty:ty, $tt:tt) => {
        impl<$tt> std::ops::Mul<MultivariatePolynomial<Self, O>> for $ty
        where
            O: MonomialOrder<Self>,
        {
            type Output = MultivariatePolynomial<Self, O>;

            fn mul(self, rhs: Self::Output) -> Self::Output {
                Self::Output::constant(None, self) * rhs
            }
        }
    };
}

impl_mul!(Real, O);
impl_mul!(Integer, O);
impl_mul!(Rational, O);

impl<R, O> std::ops::Mul<R> for MultivariatePolynomial<R, O>
where
    R: Ring,
    O: MonomialOrder<R>,
{
    type Output = Self;

    fn mul(self, rhs: R) -> Self {
        self * Self::constant(None, rhs)
    }
}
