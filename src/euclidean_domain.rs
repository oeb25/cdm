use crate::{domain::Domain, naturals::Natural};

/// <https://en.wikipedia.org/wiki/Euclidean_domain>
pub trait EuclideanDomain:
    Domain + std::ops::Rem<Output = Self> + std::ops::Div<Output = Self>
{
    fn d(&self) -> Option<Natural>;
}

#[derive(Debug)]
pub struct ExtendedEuclideanAlgorithm<D> {
    pub r: Vec<D>,
    pub s: Vec<D>,
    pub t: Vec<D>,
    pub q: Vec<D>,
}

impl<D: EuclideanDomain + PartialEq + std::fmt::Debug> ExtendedEuclideanAlgorithm<D> {
    pub fn perform(f: &D, g: &D) -> Self {
        let mut r = vec![f.clone(), g.clone()];
        let mut s = vec![D::one(), D::zero()];
        let mut t = vec![D::zero(), D::one()];
        let mut q = vec![D::zero()];

        let mut i = 1;
        while !r[i].is_zero() {
            // assert_eq!(
            //     (r[i - 1].clone() * r[i].clone()) / r[i].clone(),
            //     r[i - 1].clone()
            // );

            q.push(r[i - 1].clone() / r[i].clone());

            // println!("{:?}", (r[i - 1].clone(), r[i].clone()));

            r.push(r[i - 1].clone() - q[i].clone() * r[i].clone());
            s.push(s[i - 1].clone() - q[i].clone() * s[i].clone());
            t.push(t[i - 1].clone() - q[i].clone() * t[i].clone());

            i += 1;
        }

        Self { r, s, t, q }
    }
    pub fn gcd(&self) -> &D {
        &self.r[self.r.len() - 2]
    }
    pub fn property(f: D, g: D) -> bool {
        let res = Self::perform(&f, &g);
        let gcd = res.gcd();

        (f % gcd.clone()).is_zero() && (g % gcd.clone()).is_zero()
    }
}
