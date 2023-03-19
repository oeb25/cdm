use std::cmp::Ordering;

use itertools::Itertools;

use crate::{
    field::Field,
    identity::{Addition, Identity, Multiplication},
    multivariate_polynomials::MultivariatePolynomial,
    num_to_superscript,
    rationals::rational,
    Group, Natural, Rational, Ring,
};

pub trait MonomialOrder<F>: Clone + std::fmt::Debug {
    fn ord(&self, l: &Monomial<F, Self>, r: &Monomial<F, Self>) -> Ordering;
}

/// Pure lexicographic order.
#[derive(Debug, Clone)]
pub struct PLex(pub Vec<usize>);

impl Default for PLex {
    fn default() -> Self {
        Self((0..10).collect())
    }
}

impl<F> MonomialOrder<F> for PLex {
    fn ord(&self, l: &Monomial<F, Self>, r: &Monomial<F, Self>) -> Ordering {
        plex(Some(&self.0), l, r)
    }
}

pub fn plex<F, O>(ord: Option<&[usize]>, l: &Monomial<F, O>, r: &Monomial<F, O>) -> Ordering
where
    O: MonomialOrder<F>,
{
    ord.map(|ord| ord.to_vec())
        .unwrap_or_else(|| PLex::default().0)
        .into_iter()
        .map(|idx| {
            let l = l.powers().get(idx).cloned().unwrap_or(0);
            let r = r.powers().get(idx).cloned().unwrap_or(0);

            l.cmp(&r)
        })
        .find(|o: &Ordering| !o.is_eq())
        .unwrap_or(Ordering::Equal)
}

#[test]
fn test_plex() {
    let [x, y, z] = MultivariatePolynomial::<Rational, _>::init(PLex::default());
    let f = rational(4.) * x(1) * y(1) * z(2) + rational(4.) * x(3) - rational(5.) * y(4)
        + rational(7.) * x(1) * y(2) * z(1);

    assert_eq!(format!("{f:?}"), "4x³ + 7xy²z + 4xyz² + -5y⁴");
}

/// Graded lexicographic order.
#[derive(Debug, Clone)]
pub struct GrLex(pub Vec<usize>);

impl Default for GrLex {
    fn default() -> Self {
        Self((0..10).collect())
    }
}

impl<F> MonomialOrder<F> for GrLex {
    fn ord(&self, l: &Monomial<F, Self>, r: &Monomial<F, Self>) -> Ordering {
        grlex(Some(&self.0), l, r)
    }
}

pub fn grlex<F, O>(ord: Option<&[usize]>, l: &Monomial<F, O>, r: &Monomial<F, O>) -> Ordering
where
    O: MonomialOrder<F>,
{
    l.powers()
        .iter()
        .sum::<Natural>()
        .cmp(&r.powers().iter().sum())
        .then_with(|| plex(ord, l, r))
}

#[test]
fn test_grlex() {
    let [x, y, z] = MultivariatePolynomial::<Rational, _>::init(GrLex::default());
    let f = rational(4.) * x(1) * y(1) * z(2) + rational(4.) * x(3) - rational(5.) * y(4)
        + rational(7.) * x(1) * y(2) * z(1);

    assert_eq!(format!("{f:?}"), "7xy²z + 4xyz² + -5y⁴ + 4x³");
}

/// Graded reverse lexicographic order.
#[derive(Debug, Clone)]
pub struct TLex(pub Vec<usize>);

impl Default for TLex {
    fn default() -> Self {
        Self((0..10).collect())
    }
}

impl<F> MonomialOrder<F> for TLex
where
    F: Group,
{
    fn ord(&self, l: &Monomial<F, Self>, r: &Monomial<F, Self>) -> Ordering {
        tlex(Some(&self.0), l, r)
    }
}

pub fn tlex<F, O>(ord: Option<&[usize]>, l: &Monomial<F, O>, r: &Monomial<F, O>) -> Ordering
where
    F: Group,
    O: MonomialOrder<F>,
{
    l.powers()
        .iter()
        .sum::<Natural>()
        .cmp(&r.powers().iter().sum())
        .then_with(|| {
            ord.map(|ord| ord.to_vec())
                .unwrap_or_else(|| PLex::default().0)
                .into_iter()
                .rev()
                .map(|idx| {
                    let l = l.powers().get(idx).cloned().unwrap_or(0);
                    let r = r.powers().get(idx).cloned().unwrap_or(0);

                    r.cmp(&l)
                })
                .find(|o| !o.is_eq())
                .unwrap_or(Ordering::Equal)
        })
}

#[test]
fn test_tlex() {
    let [x, y, z] = MultivariatePolynomial::<Rational, _>::init(TLex::default());
    let f = rational(4.) * x(1) * y(1) * z(2) + rational(4.) * x(3) - rational(5.) * y(4)
        + rational(7.) * x(1) * y(2) * z(1);

    assert_eq!(format!("{f:?}"), "-5y⁴ + 7xy²z + 4xyz² + 4x³");
}
#[derive(Clone)]
pub enum Monomial<F, O: MonomialOrder<F>> {
    Constant(F),
    Mono {
        ord: O,
        coef: F,
        powers: Vec<Natural>,
    },
}
impl<F: Default, O: MonomialOrder<F>> Default for Monomial<F, O> {
    fn default() -> Self {
        Self::Constant(F::default())
    }
}

impl<F, O> PartialEq for Monomial<F, O>
where
    F: PartialEq,
    O: MonomialOrder<F>,
{
    fn eq(&self, other: &Self) -> bool {
        use itertools::EitherOrBoth;
        self.coef() == other.coef()
            && self
                .powers()
                .iter()
                .zip_longest(other.powers().iter())
                .all(|ps| match ps {
                    EitherOrBoth::Both(l, r) => l == r,
                    EitherOrBoth::Left(p) | EitherOrBoth::Right(p) => *p == 0,
                })
    }
}

impl<F, O: MonomialOrder<F>> Monomial<F, O> {
    pub fn new(ord: O, coef: impl Into<F>, powers: Vec<Natural>) -> Self {
        Self::Mono {
            ord,
            coef: coef.into(),
            powers,
        }
    }
    pub fn try_new(ord: Option<O>, coef: impl Into<F>, powers: Vec<Natural>) -> Option<Self> {
        if let Some(ord) = ord {
            return Some(Self::new(ord.clone(), coef.into(), powers));
        }
        if powers.is_empty() || powers.iter().all(|c| *c == 0) {
            Some(Self::constant(None, coef.into()))
        } else {
            None
        }
    }
    pub fn coef(&self) -> &F {
        match self {
            Self::Constant(coef) | Self::Mono { coef, .. } => coef,
        }
    }
    pub fn ord(&self) -> Option<&O> {
        match self {
            Self::Constant(_) => None,
            Self::Mono { ord, .. } => Some(ord),
        }
    }
    pub fn constant(ord: Option<O>, coef: impl Into<F>) -> Self {
        if let Some(ord) = ord {
            Self::new(ord, coef, vec![])
        } else {
            Self::Constant(coef.into())
        }
    }
    pub fn zero(ord: Option<O>) -> Self
    where
        F: Identity<Addition>,
    {
        Self::constant(ord, <F as Identity<Addition>>::identity())
    }
    pub fn is_zero(&self) -> bool
    where
        F: Identity<Addition>,
    {
        <F as Identity<Addition>>::is_identity(self.coef())
    }
    pub fn powers(&self) -> &[Natural] {
        match self {
            Self::Constant(_) => &[],
            Self::Mono { powers, .. } => powers,
        }
    }
    pub fn without_coef(&self) -> Self
    where
        F: Identity<Multiplication>,
    {
        self.map_coef(|_| F::identity())
    }
    pub fn map_coef(&self, f: impl FnOnce(&F) -> F) -> Self {
        match self {
            Self::Constant(coef) => Self::Constant(f(coef)),
            Self::Mono { ord, coef, powers } => Self::Mono {
                ord: ord.clone(),
                coef: f(coef),
                powers: powers.clone(),
            },
        }
    }
    pub fn ord_from<'a>(&'a self, rhs: &'a Self) -> Option<&'a O> {
        match (self, rhs) {
            (Self::Constant(_), Self::Constant(_)) => None,
            (Self::Mono { ord, .. }, _) => Some(ord),
            (_, Self::Mono { ord, .. }) => Some(ord),
        }
    }
    pub fn div(&self, rhs: &Self) -> Option<Self>
    where
        F: Field + std::fmt::Debug,
    {
        // eprintln!("try div {self:?}/{rhs:?}");
        use itertools::EitherOrBoth;
        let powers = self
            .powers()
            .iter()
            .zip_longest(rhs.powers())
            .map(|ps| match ps {
                EitherOrBoth::Both(l, r) => {
                    if l >= r {
                        Some(*l - *r)
                    } else {
                        None
                    }
                }
                EitherOrBoth::Left(l) => Some(*l),
                EitherOrBoth::Right(r) => {
                    if *r == 0 {
                        Some(0)
                    } else {
                        None
                    }
                }
            })
            .collect::<Option<_>>()?;

        Self::try_new(
            self.ord_from(rhs).cloned(),
            self.coef().clone() / rhs.coef().clone(),
            powers,
        )
    }
}

impl<F, O> std::fmt::Debug for Monomial<F, O>
where
    F: std::fmt::Debug + Identity<Addition> + Identity<Multiplication>,
    O: MonomialOrder<F>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if Identity::<Addition>::is_identity(self.coef()) {
            write!(f, "0")
        } else if Identity::<Multiplication>::is_identity(self.coef())
            && self.powers().iter().all(|p| *p == 0)
        {
            write!(f, "1")
        } else {
            write!(
                f,
                "{}{}",
                if Identity::<Multiplication>::is_identity(self.coef()) {
                    "".to_string()
                } else {
                    format!("{:?}", self.coef())
                },
                self.powers()
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, p)| {
                        if *p == 0 {
                            None
                        } else if Identity::<Multiplication>::is_identity(p) {
                            Some(format!("{}", ["x", "y", "z", "v", "w"][idx]))
                        } else {
                            Some(format!(
                                "{}{}",
                                ["x", "y", "z", "v", "w"][idx],
                                num_to_superscript(*p as _)
                            ))
                        }
                    })
                    .format("")
            )
        }
    }
}

impl<F, O> std::ops::Mul for Monomial<F, O>
where
    F: Ring,
    O: MonomialOrder<F>,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self::try_new(
            self.ord_from(&rhs).cloned(),
            self.coef().clone() * rhs.coef().clone(),
            self.powers()
                .iter()
                .cloned()
                .zip_longest(rhs.powers().iter().cloned())
                .map(|p| p.reduce(|l, r| l + r))
                .collect(),
        )
        .unwrap()
    }
}
impl<F, O> std::ops::Mul<F> for Monomial<F, O>
where
    F: Ring,
    O: MonomialOrder<F>,
{
    type Output = Self;

    fn mul(self, rhs: F) -> Self {
        Self::try_new(
            self.ord().cloned(),
            self.coef().clone() * rhs,
            self.powers().to_vec(),
        )
        .unwrap()
    }
}

impl<F, O> std::ops::Div<F> for Monomial<F, O>
where
    F: Field,
    O: MonomialOrder<F>,
{
    type Output = Self;

    fn div(self, rhs: F) -> Self {
        Self::try_new(
            self.ord().cloned(),
            self.coef().clone() / rhs,
            self.powers().to_vec(),
        )
        .unwrap()
    }
}
