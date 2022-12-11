use crate::{domain::Domain, naturals::Natural};

/// <https://en.wikipedia.org/wiki/Euclidean_domain>
pub trait EuclideanDomain:
    Domain + std::ops::Rem<Output = Self> + std::ops::Div<Output = Self>
{
    fn d(&self) -> Option<Natural>;
}
