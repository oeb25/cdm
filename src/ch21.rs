//! # Gröbner bases

use itertools::Itertools;
use tracing::debug;

use crate::{
    field::Field, latex::ToLatex, mono::MonomialOrder,
    multivariate_polynomials::MultivariatePolynomial, Group, Ring,
};

/// Algorithm 21.11 Multivariate division with remainder.
pub fn multivariate_division_with_remainder<F: Field, O: MonomialOrder<F>>(
    f: &MultivariatePolynomial<F, O>,
    fs: &[MultivariatePolynomial<F, O>],
) -> MultivariateDivisionWithRemainder<F, O> {
    let span = tracing::span!(
        tracing::Level::DEBUG,
        "MVDWR",
        f = format!("{:?}", f),
        fs = format!("{:?}", fs)
    );
    let _enter = span.enter();

    let mut result = MultivariateDivisionWithRemainder::new(f.clone(), fs.to_vec());

    let mut r = MultivariatePolynomial::<F, O>::zero();
    let mut p = f.clone();
    let mut q = fs
        .iter()
        .map(|_| MultivariatePolynomial::<F, O>::zero())
        .collect_vec();

    // eprintln!("{:?}", p.normalized());

    result.add_row(p.clone(), q.to_vec(), r.clone());

    while !p.is_zero() {
        debug!("in here... {p:?}");

        let lt_p = p.leading_term();

        let mut found = false;
        for (i, fi) in fs.iter().enumerate() {
            let lt_f = fi.leading_term();

            debug!("p = {p:?}, lt(p) = {lt_p:?}");
            debug!("lt(f) = {lt_f:?}");

            if let Some(ratio) = lt_p.div(&lt_f) {
                found = true;

                debug!("RATIO! {ratio:?}");

                let ratio: MultivariatePolynomial<_, O> = ratio.into();

                q[i] = q[i].clone() + ratio.clone();

                debug!(
                    "subtracting ({ratio:?})*({fi:?}) = {:?}",
                    ratio.clone() * fi.clone()
                );
                p = p.clone() - ratio * fi.clone();

                break;
            }
        }

        debug!("found = {found:?}");

        if !found {
            r = r + lt_p.clone().into();
            debug!("Reducing p = {p:?} by lt(p) = {lt_p:?}");
            p = p - lt_p.into();
            debug!("p is thus {p:?}");
        }

        result.add_row(p.clone(), q.to_vec(), r.clone());
    }

    // eprintln!("{}", result.ascii_table(ord));

    result
}

pub struct MultivariateDivisionWithRemainder<F, O: MonomialOrder<F>> {
    pub f: MultivariatePolynomial<F, O>,
    pub fs: Vec<MultivariatePolynomial<F, O>>,
    pub rows: Vec<Vec<MultivariatePolynomial<F, O>>>,
}

impl<F: Field, O: MonomialOrder<F>> MultivariateDivisionWithRemainder<F, O> {
    pub fn new(f: MultivariatePolynomial<F, O>, fs: Vec<MultivariatePolynomial<F, O>>) -> Self {
        Self {
            f,
            fs,
            rows: vec![],
        }
    }

    fn add_row(
        &mut self,
        p: MultivariatePolynomial<F, O>,
        mut q: Vec<MultivariatePolynomial<F, O>>,
        r: MultivariatePolynomial<F, O>,
    ) {
        q.insert(0, p);
        q.push(r);
        self.rows.push(q);
    }

    pub fn latex_table(&self) -> String
    where
        F: ToLatex + Field,
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
            write!(buf, "{}", e[0].to_latex());
            for p in &e[1..] {
                write!(buf, "& {}", p.to_latex());
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

    pub fn ascii_table(&self) -> comfy_table::Table {
        use comfy_table::{
            modifiers::UTF8_ROUND_CORNERS, presets::UTF8_FULL, ContentArrangement, Table,
        };

        let mut table = Table::new();
        table
            .load_preset(UTF8_FULL)
            .apply_modifier(UTF8_ROUND_CORNERS)
            .set_content_arrangement(ContentArrangement::Dynamic);

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
                std::iter::once(format!("{:?}", p))
                    .chain(q.iter().map(|q| format!("{:?}", q)))
                    .chain(std::iter::once(format!("{:?}", r))),
            );
        }

        table
    }

    pub(crate) fn remainder(&self) -> MultivariatePolynomial<F, O> {
        self.rows.last().unwrap().last().unwrap().clone()
    }
}

#[derive(Debug)]
pub struct GroebnerBasis<F: Ring, O: MonomialOrder<F>> {
    gs: Vec<MultivariatePolynomial<F, O>>,
}

/// Algorithm 21.33 Gröbner basis computation.
pub fn buchbergers_algorithm<F, O>(fs: &[MultivariatePolynomial<F, O>]) -> GroebnerBasis<F, O>
where
    F: Field,
    O: MonomialOrder<F>,
{
    let mut basis = fs.to_vec();
    basis.dedup();

    loop {
        let mut new = vec![];

        for (i, f1) in basis.iter().enumerate() {
            for f2 in &basis[i + 1..] {
                let s = f1.s_polynomial(f2);
                let result = multivariate_division_with_remainder(&s, &basis);

                let r = result.remainder();

                if !r.is_zero() && !new.contains(&r) && !basis.contains(&r) {
                    if !multivariate_division_with_remainder(&r, &new)
                        .remainder()
                        .is_zero()
                    {
                        new.push(r);
                    }
                }
            }
        }

        if new.is_empty() {
            break;
        }

        basis.extend_from_slice(&new);
    }

    GroebnerBasis { gs: basis }
}

/// Definition 21.37
pub fn minimize_groebner_basis<F, O>(basis: &mut GroebnerBasis<F, O>)
where
    F: Field,
    O: MonomialOrder<F>,
{
    let mut i = 0;

    while i < basis.gs.len() {
        basis.gs.swap(0, i);

        let f = basis.gs[0].leading_term().into();
        let rest = basis.gs[1..]
            .iter()
            .map(|g| g.leading_term().into())
            .collect_vec();

        let r = multivariate_division_with_remainder(&f, &rest).remainder();

        if r.is_zero() {
            basis.gs.remove(0);
        } else {
            i += 1;
        }
    }

    for g in &mut basis.gs {
        *g = g.minimize();
    }
}

#[cfg(test)]
mod tests {
    use crate::{mono::PLex, multivariate_polynomials::MultivariatePolynomial, Rational, Ring};

    use super::multivariate_division_with_remainder;

    #[test]
    fn multi_poly_test_1() {
        let [x, y] = MultivariatePolynomial::init(PLex(vec![0, 1]));

        let f = x(2) * y(1) + x(1) * y(2) + y(2);
        let f1 = x(1) * y(1) - MultivariatePolynomial::one();
        let f2 = y(2) - MultivariatePolynomial::one();

        multivariate_division_with_remainder::<Rational, _>(&f, &[f1, f2]);
    }

    #[test]
    fn multi_poly_test_2() {
        tracing_subscriber::fmt::fmt()
            .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
            .without_time()
            .init();

        let [x, y] = MultivariatePolynomial::<Rational, PLex>::init(PLex(vec![0, 1]));

        let f = x(2) * y(1) + y(10) + MultivariatePolynomial::one();
        let f1 = y(1) + MultivariatePolynomial::one();
        let f2 = y(7) - MultivariatePolynomial::one();
        let f3 = x(1) * y(5) + MultivariatePolynomial::constant(None, Rational::approximate(5.0));

        multivariate_division_with_remainder(&f, &[f1, f2, f3]);
    }
}
