use std::cmp::Ordering;

use itertools::Itertools;

use crate::{
    field::Field,
    identity::{Addition, Identity, Multiplication},
    num_to_superscript, Natural, Ring,
};

pub trait MonomialOrder<F> {
    fn ord(&self, l: &Monomial<F>, r: &Monomial<F>) -> Ordering;
}

pub struct PLex(pub Vec<usize>);

impl Default for PLex {
    fn default() -> Self {
        Self((0..10).collect())
    }
}

impl<F> MonomialOrder<F> for PLex {
    fn ord(&self, l: &Monomial<F>, r: &Monomial<F>) -> Ordering {
        self.0
            .iter()
            .cloned()
            .map(|idx| {
                let l = l.powers.get(idx).cloned().unwrap_or(0);
                let r = r.powers.get(idx).cloned().unwrap_or(0);

                r.cmp(&l)
            })
            .find(|o: &Ordering| !o.is_eq())
            .unwrap_or(Ordering::Equal)
    }
}

#[derive(Clone, Default)]
pub struct Monomial<F> {
    coef: F,
    powers: Vec<Natural>,
}

impl<F> PartialEq for Monomial<F>
where
    F: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        // TODO
        self.coef == other.coef
            && self
                .powers
                .iter()
                .zip_longest(other.powers.iter())
                .map(|ps| match ps {
                    itertools::EitherOrBoth::Both(l, r) => l == r,
                    itertools::EitherOrBoth::Left(p) | itertools::EitherOrBoth::Right(p) => *p == 0,
                })
                .all(|x| x)
    }
}

impl<F> Monomial<F> {
    pub fn new(coef: impl Into<F>, powers: Vec<Natural>) -> Self {
        Monomial {
            coef: coef.into(),
            powers: powers.into_iter().map_into().collect(),
        }
    }
    pub fn constant(coef: impl Into<F>) -> Self {
        Monomial {
            coef: coef.into(),
            powers: vec![],
        }
    }
    pub fn zero() -> Self
    where
        F: Identity<Addition>,
    {
        Monomial {
            coef: <F as Identity<Addition>>::identity(),
            powers: vec![],
        }
    }
    pub fn is_zero(&self) -> bool
    where
        F: Identity<Addition>,
    {
        self.coef == <F as Identity<Addition>>::identity()
    }
    pub fn coef(&self) -> &F {
        &self.coef
    }
    pub fn powers(&self) -> &[Natural] {
        &self.powers
    }
    pub fn without_coef(&self) -> Monomial<F>
    where
        F: Identity<Multiplication>,
    {
        Monomial {
            coef: F::identity(),
            powers: self.powers.clone(),
        }
    }
    pub fn map_coef<T>(&self, f: impl FnOnce(&F) -> T) -> Monomial<T> {
        Monomial {
            coef: f(&self.coef),
            powers: self.powers.clone(),
        }
    }
    pub fn div(&self, rhs: &Self) -> Option<Self>
    where
        F: Field + std::fmt::Debug,
    {
        // eprintln!("try div {self:?}/{rhs:?}");

        let powers = self
            .powers
            .iter()
            .zip_longest(rhs.powers.iter())
            .map(|ps| match ps {
                itertools::EitherOrBoth::Both(l, r) => {
                    if l >= r {
                        Some(*l - *r)
                    } else {
                        None
                    }
                }
                itertools::EitherOrBoth::Left(l) => Some(*l),
                itertools::EitherOrBoth::Right(r) => {
                    if *r == 0 {
                        Some(0)
                    } else {
                        None
                    }
                }
            })
            .collect::<Option<_>>();

        let res = Monomial::new(self.coef.clone() / rhs.coef.clone(), powers.clone()?);

        // eprintln!("gives ==> ({:?}) {:?}", res, powers);

        Some(res)
    }
}

impl<F> std::fmt::Debug for Monomial<F>
where
    F: std::fmt::Debug + Identity<Addition> + Identity<Multiplication>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if Identity::<Addition>::is_identity(&self.coef) {
            write!(f, "0")
        } else if Identity::<Multiplication>::is_identity(&self.coef)
            && self.powers.iter().all(|p| *p == 0)
        {
            write!(f, "1")
        } else {
            write!(
                f,
                "{}{}",
                if Identity::<Multiplication>::is_identity(&self.coef) {
                    "".to_string()
                } else {
                    format!("{:?}", self.coef)
                },
                self.powers
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

impl<F> std::ops::Mul for Monomial<F>
where
    F: Ring,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let l = format!("{:?}", self);
        let r = format!("{:?}", rhs);

        let res = Monomial {
            coef: self.coef * rhs.coef,
            powers: self
                .powers
                .iter()
                .cloned()
                .zip_longest(rhs.powers.iter().cloned())
                .map(|p| p.reduce(|l, r| l + r))
                .collect(),
        };

        // eprintln!("MONO: {l}*{r} = {res:?}");

        res
    }
}
impl<F> std::ops::Mul<F> for Monomial<F>
where
    F: Ring,
{
    type Output = Self;

    fn mul(self, rhs: F) -> Self::Output {
        Monomial {
            coef: self.coef * rhs,
            powers: self.powers.clone(),
        }
    }
}

impl<F> std::ops::Div<F> for Monomial<F>
where
    F: Field,
{
    type Output = Self;

    fn div(self, rhs: F) -> Self::Output {
        Monomial {
            coef: self.coef / rhs,
            powers: self.powers,
        }
    }
}
