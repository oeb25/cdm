use itertools::Itertools;

use crate::{
    field::Field,
    identity::{Addition, Identity, Multiplication},
    latex::ToLatex,
    mono::{Monomial, MonomialOrder, PLex},
    Group, Natural, Rational, Ring,
};

#[derive(Clone)]
pub struct MultivariatePolynomial<F> {
    pub terms: Vec<Monomial<F>>,
}

impl<F> PartialEq for MultivariatePolynomial<F>
where
    F: Identity<Addition> + Clone + PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        // TODO
        let ord = PLex::default();

        self.terms_sorted(&ord)
            .zip_longest(other.terms_sorted(&ord))
            .map(|ps| match ps {
                itertools::EitherOrBoth::Both(l, r) => l == r,
                itertools::EitherOrBoth::Left(_) | itertools::EitherOrBoth::Right(_) => false,
            })
            .all(|x| x)
    }
}

impl<F> std::fmt::Debug for MultivariatePolynomial<F>
where
    F: std::fmt::Debug + Identity<Addition> + Identity<Multiplication>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = format!(
            "{:?}",
            self.terms.iter().filter(|t| !t.is_zero()).format(" + ")
        );
        if s.is_empty() {
            write!(f, "0")
        } else {
            write!(f, "{}", s)
        }
    }
}

impl<F> MultivariatePolynomial<F>
where
    F: Identity<Addition> + Clone,
{
    pub fn sort_by(&mut self, ord: &impl MonomialOrder<F>) {
        self.terms.sort_by(|a, b| ord.ord(a, b));
    }
    pub fn sorted_by(&self, ord: &impl MonomialOrder<F>) -> Self {
        let mut new = self.clone();
        new.sort_by(ord);
        new
    }
    pub fn leading_term(&self, ord: &impl MonomialOrder<F>) -> Monomial<F> {
        self.terms
            .iter()
            .min_by(|l, r| ord.ord(l, r))
            .cloned()
            .unwrap_or_else(|| Monomial::zero())
    }
    pub fn leading_coef(&self, ord: &impl MonomialOrder<F>) -> F {
        self.leading_term(ord).coef().clone()
    }
    pub fn leading_monomial(&self, ord: &impl MonomialOrder<F>) -> Monomial<F>
    where
        F: Identity<Multiplication>,
    {
        self.leading_term(ord).without_coef()
    }
    pub fn multi_deg(&self, ord: &impl MonomialOrder<F>) -> Vec<Natural> {
        todo!()
    }
    pub fn terms_sorted(&self, ord: &impl MonomialOrder<F>) -> impl Iterator<Item = &Monomial<F>> {
        self.terms
            .iter()
            .filter(|t| !t.is_zero())
            .sorted_by(|l, r| ord.ord(l, r))
    }
    pub fn normalized(&self) -> Self
    where
        F: Ring + Clone + std::fmt::Debug,
    {
        let ord = PLex::default();

        let terms = self
            .terms_sorted(&ord)
            .group_by(|t| t.without_coef())
            .into_iter()
            .map(|(a, b)| a * b.map(|t| t.coef()).cloned().reduce(|a, b| a + b).unwrap())
            .collect();

        MultivariatePolynomial { terms }
    }
    pub fn minimize(&self, ord: &impl MonomialOrder<F>) -> Self
    where
        F: Field,
    {
        let terms = self
            .terms_sorted(ord)
            .map(|t| t.clone() / self.leading_coef(ord))
            .collect();

        MultivariatePolynomial { terms }
    }
    pub fn s_polynomial(&self, other: &Self, ord: &impl MonomialOrder<F>) -> Self
    where
        F: Field + std::fmt::Debug,
    {
        let lcm = self
            .leading_term(ord)
            .powers()
            .iter()
            .copied()
            .zip_longest(other.leading_term(ord).powers().iter().copied())
            .map(|ps| ps.reduce(|a, b| a.max(b)))
            .collect_vec();

        let lcm = Monomial::new(F::one(), lcm.clone());

        let l = lcm
            .div(&self.leading_term(ord))
            .unwrap_or_else(|| panic!("{lcm:?}/{:?} failed", self.leading_term(ord)));
        let r = lcm
            .div(&other.leading_term(ord))
            .unwrap_or_else(|| panic!("{lcm:?}/{:?} failed", other.leading_term(ord)));

        let l = MultivariatePolynomial::from(l);
        let r = MultivariatePolynomial::from(r);

        let res = l * self.clone() - r * other.clone();

        res.normalized().sorted_by(ord)
    }
    pub fn div_mono(&self, m: &Monomial<F>) -> Option<Self>
    where
        F: Field + std::fmt::Debug,
    {
        let terms = self.terms.iter().map(|t| t.div(m)).collect::<Option<_>>()?;

        Some(MultivariatePolynomial { terms })
    }
}

impl<F> From<Monomial<F>> for MultivariatePolynomial<F> {
    fn from(value: Monomial<F>) -> Self {
        MultivariatePolynomial { terms: vec![value] }
    }
}

impl<F> std::ops::Add for MultivariatePolynomial<F>
where
    F: Ring,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        MultivariatePolynomial {
            terms: self
                .terms
                .into_iter()
                .chain(rhs.terms.into_iter())
                .collect(),
        }
        .normalized()
    }
}
impl<F> std::ops::Neg for MultivariatePolynomial<F>
where
    F: Ring + Clone,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        MultivariatePolynomial {
            terms: self
                .terms
                .into_iter()
                .map(|t| t.map_coef(|c| -c.clone()))
                .collect(),
        }
        .normalized()
    }
}
impl<F> std::ops::Sub for MultivariatePolynomial<F>
where
    F: Ring + Clone,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<F: Group> Identity<Addition> for MultivariatePolynomial<F> {
    fn identity() -> Self {
        MultivariatePolynomial { terms: vec![] }
    }
}
impl<F: Ring> Group for MultivariatePolynomial<F> {}

impl<F: Ring> std::ops::Mul for MultivariatePolynomial<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let terms = self
            .terms
            .iter()
            .cartesian_product(rhs.terms.iter())
            .map(|(a, b)| a.clone() * b.clone())
            .collect();

        MultivariatePolynomial { terms }.normalized()
    }
}

#[test]
fn multi_poly_test_1() {
    let f = MultivariatePolynomial {
        terms: vec![
            Monomial::new(1, vec![2, 1]),
            Monomial::new(1, vec![1, 2]),
            Monomial::new(1, vec![0, 2]),
        ],
    };
    let f1 = MultivariatePolynomial {
        terms: vec![Monomial::new(1, vec![1, 1]), Monomial::constant(-1)],
    };
    let f2 = MultivariatePolynomial {
        terms: vec![Monomial::new(1, vec![0, 2]), Monomial::constant(-1)],
    };

    multivariate_division_with_remainder::<Rational>(&PLex(vec![0, 1]), &f, &[f1, f2]);
}

#[test]
fn multi_poly_test_2() {
    let f = MultivariatePolynomial {
        terms: vec![
            Monomial::new(1, vec![2, 1]),
            Monomial::new(1, vec![0, 10]),
            Monomial::constant(1),
        ],
    }
    .normalized();
    let f1 = MultivariatePolynomial {
        terms: vec![Monomial::new(1, vec![1]), Monomial::constant(1)],
    }
    .normalized();
    let f2 = MultivariatePolynomial {
        terms: vec![Monomial::new(1, vec![0, 7]), Monomial::constant(-1)],
    }
    .normalized();
    let f3 = MultivariatePolynomial {
        terms: vec![Monomial::new(1, vec![1, 5]), Monomial::constant(5)],
    }
    .normalized();

    multivariate_division_with_remainder::<Rational>(&PLex(vec![0, 1]), &f, &[f1, f2, f3]);
}

pub struct MultivariateDivisionWithRemainder<F> {
    pub f: MultivariatePolynomial<F>,
    pub fs: Vec<MultivariatePolynomial<F>>,
    pub rows: Vec<Vec<MultivariatePolynomial<F>>>,
}

impl<F: Field + std::fmt::Debug> MultivariateDivisionWithRemainder<F> {
    fn new(f: MultivariatePolynomial<F>, fs: Vec<MultivariatePolynomial<F>>) -> Self {
        Self {
            f,
            fs,
            rows: vec![],
        }
    }

    fn add_row(
        &mut self,
        p: MultivariatePolynomial<F>,
        mut q: Vec<MultivariatePolynomial<F>>,
        r: MultivariatePolynomial<F>,
    ) {
        q.insert(0, p);
        q.push(r);
        self.rows.push(q);
    }

    pub fn latex_table(&self, ord: &impl MonomialOrder<F>) -> String
    where
        F: ToLatex + Identity<Addition> + Identity<Multiplication>,
    {
        use std::fmt::Write;

        let mut buf = String::new();

        write!(buf, "\\begin{{array}}{{r");
        for _ in &self.fs {
            write!(buf, "|r");
        }
        writeln!(buf, "|l}}");

        write!(buf, "\\multicolumn{{1}}{{c}}{{p}}");
        for (i, _) in self.fs.iter().enumerate() {
            write!(buf, "& \\multicolumn{{1}}{{|c}}{{q_{i}}}");
        }
        write!(buf, "& \\multicolumn{{1}}{{|c}}{{r}}");
        writeln!(buf, "\\\\ \\hline");

        for e in &self.rows {
            write!(buf, "{}", e[0].sorted_by(ord).to_latex());
            for p in &e[1..] {
                write!(buf, "& {}", p.sorted_by(ord).to_latex());
            }
            writeln!(buf, " \\\\ \\hline");
        }

        writeln!(buf, "\\end{{array}}");

        buf
        // $$
        // \begin{array}{r|r|r|r|l}
        // \multicolumn{1}{c|}{p} & \multicolumn{1}{c|}{q_1} & \multicolumn{1}{c|}{q_2} & \multicolumn{1}{c|}{q_3} & \multicolumn{1}{c}{r} \\\hline
        // - 2y^2z + 2z -2x^2y + 2xy & 0 & 0 & 0 & 0 \\\hline
        // 2z -3x^2 y +4 x y -y & 0 & 0 & y & 0 \\\hline
        // -3x^2 y +4 x y -y & 0 & 0 & y &  2z \\\hline
        // 4 x y + 3y^3 -4 y & -3 y & 0 & y &  2z \\\hline
        // 3y^3 - 4y & -3 y & 0 & y & 2z + 4 x y \\\hline
        // -4y & -3 y & 0 & y & 2z +4 x y + 3 y^3 \\\hline
        // 0 & -3 y & 0 & y & 2z +4 x y + 3 y^3  -4 y \\\hline
        // \end{array}
        // $$
    }

    pub fn ascii_table(&self, ord: &impl MonomialOrder<F>) -> comfy_table::Table {
        let mut table = comfy_table::Table::new();
        table.set_header(
            std::iter::once("p".to_string())
                // .chain(fs.iter().enumerate().map(|(i, _)| format!("q{i}")))
                .chain(self.fs.iter().enumerate().map(|(i, f)| format!("{f:?}")))
                .chain(std::iter::once("r".to_string())),
        );

        for row in &self.rows {
            let p = &row[0];
            let q = &row[1..row.len() - 1];
            let r = row.last().unwrap();

            table.add_row(
                std::iter::once(format!("{:?}", p.sorted_by(ord)))
                    .chain(q.iter().map(|q| format!("{:?}", q.sorted_by(ord))))
                    .chain(std::iter::once(format!("{:?}", r.sorted_by(ord)))),
            );
        }

        table
    }

    pub(crate) fn remainder(&self) -> MultivariatePolynomial<F> {
        self.rows.last().unwrap().last().unwrap().clone()
    }
}

/// Algorithm 21.11
pub fn multivariate_division_with_remainder<F: Field + std::fmt::Debug>(
    ord: &impl MonomialOrder<F>,
    f: &MultivariatePolynomial<F>,
    fs: &[MultivariatePolynomial<F>],
) -> MultivariateDivisionWithRemainder<F> {
    let mut result = MultivariateDivisionWithRemainder::new(f.clone(), fs.to_vec());

    let mut r = MultivariatePolynomial::<F>::zero();
    let mut p = f.clone();
    let mut q = fs
        .iter()
        .map(|_| MultivariatePolynomial::<F>::zero())
        .collect_vec();

    // eprintln!("{:?}", p.normalized());

    result.add_row(p.clone(), q.to_vec(), r.clone());

    while !p.is_zero() {
        let lt_p = p.leading_term(ord);

        let mut found = false;
        for (i, fi) in fs.iter().enumerate() {
            let lt_f = fi.leading_term(ord);

            if let Some(ratio) = lt_p.div(&lt_f) {
                found = true;

                // eprintln!("RATIO! {ratio:?}");

                let ratio: MultivariatePolynomial<_> = ratio.into();

                q[i] = q[i].clone() + ratio.clone();

                // eprintln!(
                //     "subtracting ({ratio:?})*({fi:?}) = {:?}",
                //     ratio.clone() * fi.clone()
                // );
                p = p.clone() - ratio * fi.clone();

                break;
            }
        }

        if !found {
            r = r + lt_p.clone().into();
            p = p - lt_p.into();
        }

        p = p.normalized();

        result.add_row(p.clone(), q.to_vec(), r.clone());
    }

    // eprintln!("{}", result.ascii_table(ord));

    result
}
