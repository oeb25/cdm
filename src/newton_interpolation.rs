use std::{collections::HashSet, fmt::Debug, hash::Hash};

use crate::{field::Field, polynomials::Polynomial, Group};

pub struct NewtonInterpolation<'a, F: Field> {
    samples: &'a [(F, F)],
    stuff: Vec<Vec<Option<F>>>,
}

#[derive(Debug)]
pub enum NewtonInterpolationError {
    SamplesNotUnique,
}

impl<'a, F> NewtonInterpolation<'a, F>
where
    F: Field + Hash + Eq + Debug,
{
    fn u(&self, i: usize) -> F {
        self.samples[i].0.clone()
    }
    fn v(&self, i: usize) -> F {
        self.samples[i].1.clone()
    }
    fn set_f(&mut self, i: usize, j: usize, e: F) {
        self.stuff[i][j] = Some(e);
    }
    fn f(&mut self, i: usize, j: usize) -> F {
        if i == 0 {
            return self.v(j);
        }
        if let Some(e) = self.stuff[i][j].clone() {
            e
        } else {
            // let e = (self.f(i - 1, j) - self.f(i - 1, i - 1)) / (self.u(j) - self.u(i - 1));
            let e = (self.f(i - 1, j) - self.f(i - 1, i - 1)) / (self.u(j) - self.u(i - 1));
            self.stuff[i][j] = Some(e.clone());
            e
        }
    }
    pub fn run(samples: &'a [(F, F)]) -> Result<Polynomial<F>, NewtonInterpolationError> {
        let n = samples.len();

        if n == 0 {
            return Ok(Polynomial::zero());
        }

        let us: HashSet<_> = samples.iter().map(|(u, _)| u).collect();

        if samples.len() != us.len() {
            return Err(NewtonInterpolationError::SamplesNotUnique);
        }

        let stuff: Vec<Vec<Option<F>>> = vec![vec![None; n]; n];

        let mut ni = NewtonInterpolation { samples, stuff };

        for (i, (_u, v)) in samples.iter().enumerate() {
            ni.stuff[0][i] = Some(v.clone());
        }

        for i in 1..n {
            for j in i..n {
                assert!(ni.stuff[i][j].is_none());
                println!(
                    "f_{i}({j}) = (f_{is}({j}) - f_{is}({is})) / (u_{j} - u_{is}) = ({:?} - {:?}) / ({:?} - {:?}) = {:?}",
                    ni.stuff[i - 1][j].clone().unwrap(), ni.stuff[i - 1][i - 1].clone().unwrap()
                    , samples[j].0.clone(), samples[i - 1].0.clone(),
                    (ni.stuff[i - 1][j].clone().unwrap() - ni.stuff[i - 1][i - 1].clone().unwrap())
                        / (samples[j].0.clone() - samples[i - 1].0.clone()),
                    is = i - 1,

                );

                ni.stuff[i][j] = Some(
                    (ni.stuff[i - 1][j].clone().unwrap() - ni.stuff[i - 1][i - 1].clone().unwrap())
                        / (samples[j].0.clone() - samples[i - 1].0.clone()),
                );
            }
        }
        // n(n-1)/2 * 3

        println!("{:?}", ni.stuff);

        // let mut p = Polynomial::new(vec![ni.f(n - 1, n - 1)]);
        let mut p = Polynomial::new(vec![ni.stuff[n - 1][n - 1].clone().unwrap()]);
        println!("f := {p:?}");
        for i in (0..n - 1).rev() {
            // p = p.times_x(1) - p.scale(&ni.u(i)) + Polynomial::new(vec![ni.f(i, i)]);
            let p1 = p.times_x(1) - p.scale(&ni.u(i))
                + Polynomial::new(vec![ni.stuff[i][i].clone().unwrap()]);
            println!(
                "f := (x - {:?})({:?}) + {:?} = {p1:?}",
                ni.u(i),
                p,
                &ni.stuff[i][i].clone().unwrap()
            );
            p = p1;
            // 2i
        }
        // sum(i=0..n) 2i = n(n-1)/2 * 2 = n(n-1)
        // n*(n-1)

        // n(n-1)/2 * 3 + n*(n-1)

        Ok(p)
    }
}

fn newton<F, Map>(samples: &[F], map: Map) -> Polynomial<F>
where
    F: Field + Debug,
    Map: Fn(usize) -> F,
{
    if samples.is_empty() {
        return Polynomial::zero();
    }

    // (x - u0)g = x*g - u0*g
    let g = newton(
        &samples[1..],
        Box::new(|i: usize| (map(i + 1) - map(0)) / (samples[i + 1].clone() - samples[0].clone()))
            as Box<dyn Fn(usize) -> F>,
    );

    g.times_x(1) - g.scale(&samples[0]) + Polynomial::new(vec![map(0)])
}

fn newton_slice<F, Map>(samples: &[F], offset: usize, map: Map) -> Polynomial<F>
where
    F: Field + Debug,
    Map: Fn(usize) -> F,
{
    if samples.len() == offset {
        return Polynomial::zero();
    }

    // (x - u0)g = x*g - u0*g
    let g = newton_slice(
        samples,
        offset + 1,
        Box::new(|i: usize| {
            (map(i + 1) - map(0)) / (samples[offset + i + 1].clone() - samples[offset].clone())
        }) as Box<dyn Fn(usize) -> F>,
    );

    g.times_x(1) - g.scale(&samples[offset]) + Polynomial::new(vec![map(0)])
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use proptest::{prelude::prop, prop_assume, prop_compose, proptest};

    use crate::{
        count_ops::{self, CountOps},
        Finite, Rational,
    };

    use super::{NewtonInterpolation, NewtonInterpolationError};

    #[test]
    fn example_a() -> Result<(), NewtonInterpolationError> {
        // let samples: [(CountOps<Rational>, CountOps<Rational>); 3] = [(0, 1), (1, 2), (2, 4)]
        let samples: [(Finite<5>, Finite<5>); 3] = [(0, 1), (1, 2), (2, 4)].map(|(u, v)| {
            (
                Finite::from(u as i128).into(),
                Finite::from(v as i128).into(),
            )
        });

        let f = NewtonInterpolation::run(&samples)?;
        // let f = newton_slice(&samples.map(|(u, v)| u), 0, |i| samples[i].1);

        println!("{f:?}");

        for (u, v) in samples {
            assert_eq!(f.evaluate_at(u), v);
        }

        // panic!();

        Ok(())
    }
    #[test]
    fn example_b() -> Result<(), NewtonInterpolationError> {
        let samples: [(CountOps<Finite<7>>, CountOps<Finite<7>>); 3] = [
            (0, 1),
            (1, 5),
            (6, 2),
            // (10, 3),
            // (9, 4),
            // (11, 1),
            // (12, 5),
            // (17, 2),
            // (21, 3),
            // (14, 4),
        ]
        .map(|(u, v)| {
            (
                Finite::from(u as i128).into(),
                Finite::from(v as i128).into(),
            )
        });

        count_ops::reset();
        let f = NewtonInterpolation::run(&samples)?;
        // let f = newton(&samples.map(|(u, v)| u), |i| samples[i].1);
        let ops = count_ops::get_counts();
        // assert!((ops.additions + ops.multiplications) as usize <= (5 * samples.len().pow(2)) / 2);

        for (u, v) in samples {
            assert_eq!(f.evaluate_at(u), v, "f({u:?}) != {v:?}");
        }

        Ok(())
    }
    #[test]
    fn example_c() -> Result<(), NewtonInterpolationError> {
        // let samples: [(CountOps<Rational>, CountOps<Rational>); 3] = [(0, 1), (1, 2), (2, 4)]
        let samples: [(Rational, Rational); 4] = [(-5, -2), (-1, 6), (0, -1), (2, 3)]
            .map(|(u, v)| (Rational::from(u).into(), Rational::from(v).into()));

        let f = NewtonInterpolation::run(&samples)?;
        println!("{f:?}");

        for (u, v) in samples {
            assert_eq!(f.evaluate_at(u), v);
        }

        Ok(())
    }

    prop_compose! {
        fn samples()(n in 0..101usize)(us in prop::collection::vec(1..10u128, n), vs in prop::collection::vec(0..10u128, n)) -> Vec<(u128, u128)> {
            us.into_iter().scan(0, |state,x| {
                *state = *state + x;
                Some(*state)
            }).zip(vs.into_iter()).collect_vec()
        }
    }

    proptest! {
        #[test]
        fn newton_interp(uv in samples()) {
            let uv: Vec<(CountOps<Finite<101>>, CountOps<Finite<101>>)> = uv.into_iter().map(|(u, v)| (Finite::from(u).into(), Finite::from(v).into())).collect();

            count_ops::reset();
            let f = NewtonInterpolation::run(&uv);
            prop_assume!(f.is_ok());
            let f = f.unwrap();
            // let f = newton_slice(&uv.iter().map(|(u, v)| u.clone()).collect_vec(), 0, |i| uv[i].1);
            let ops = count_ops::get_counts();
            // assert!((ops.additions + ops.multiplications) as usize <= (5 * uv.len().pow(2)) / 2);

            for (u, v) in uv {
                assert_eq!(f.evaluate_at(u), v, "f({u:?}) != {v:?}");
            }
        }
    }
}
