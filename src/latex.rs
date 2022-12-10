use crate::{
    identity::{Addition, Identity, Multiplication},
    mono::Monomial,
    multivariate_polynomials::MultivariatePolynomial,
    Finite, Natural, Rational,
};
use itertools::Itertools;

pub trait ToLatex {
    fn to_latex(&self) -> String;
}

impl<F> ToLatex for MultivariatePolynomial<F>
where
    F: ToLatex + Identity<Addition> + Identity<Multiplication>,
{
    fn to_latex(&self) -> String {
        let s = format!(
            "{}",
            self.terms
                .iter()
                .filter(|t| !t.is_zero())
                .map(|t| t.to_latex())
                .format(" + ")
        );
        if s.is_empty() {
            format!("0")
        } else {
            format!("{}", s)
        }
    }
}

impl<F> ToLatex for Monomial<F>
where
    F: ToLatex + Identity<Addition> + Identity<Multiplication>,
{
    fn to_latex(&self) -> String {
        use std::fmt::Write;

        let mut f = String::new();

        if Identity::<Addition>::is_identity(self.coef()) {
            write!(f, "0");
        } else if Identity::<Multiplication>::is_identity(self.coef())
            && self.powers().iter().all(|p| *p == 0)
        {
            write!(f, "1");
        } else {
            write!(
                f,
                "{}{}",
                if Identity::<Multiplication>::is_identity(self.coef()) {
                    "".to_string()
                } else {
                    self.coef().to_latex()
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
                                "{}^{{{:?}}}",
                                ["x", "y", "z", "v", "w"][idx],
                                Natural::from(*p)
                            ))
                        }
                    })
                    .format("")
            );
        }
        f
    }
}

impl ToLatex for Rational {
    fn to_latex(&self) -> String {
        if self.denom == 1 {
            format!("{:?}", self.num)
        } else if self.denom == 1 {
            format!("-{:?}", self.num)
        } else if self.num == self.num.abs() {
            format!("\\nicefrac{{{:?}}}{{{:?}}}", self.num.abs(), self.denom)
        } else {
            format!("-\\nicefrac{{{:?}}}{{{:?}}}", self.num.abs(), self.denom)
        }
    }
}

impl<const N: u128> ToLatex for Finite<N> {
    fn to_latex(&self) -> String {
        format!("{self:?}")
    }
}
