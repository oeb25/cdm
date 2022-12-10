use crate::{polynomials::Polynomial, ring::Ring};

/// A commutative rings without zero-divisors is called an integral domain (also just called a domain).
pub trait Domain: Ring {}

impl<F> Domain for Polynomial<F> where F: Ring + Clone {}
