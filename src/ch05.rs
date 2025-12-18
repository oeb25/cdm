//! # Modular algorithms and interpolation

use tracing::debug;

use crate::{ch03::extended_euclidean_algorithm, euclidean_domain::EuclideanDomain};

/// Algorithm 5.4 Chinese Remainder Algorithm (CRA).
pub fn chinese_remainder_algorithm<R: EuclideanDomain + PartialOrd>(ms: &[R], v: &[R]) -> R {
    let scope = tracing::span!(
        tracing::Level::DEBUG,
        "(5.4) CRA",
        ms = format!("{ms:?}"),
        v = format!("{v:?}"),
    );
    let _enter = scope.enter();

    let mut c = vec![];
    let m = ms.iter().cloned().reduce(|m0, m1| m0 * m1).unwrap();
    debug!("m = {m:?}");
    for (i, (mi, vi)) in ms.iter().zip(v).enumerate() {
        let res = extended_euclidean_algorithm(&(m.clone() / mi.clone()), mi);
        debug!("EEA on {:?} and {:?}", m.clone() / mi.clone(), mi);
        println!("{}, {:?}", res, &res.s[res.s.len() - 2]);
        let ci = (vi.clone() * res.s[res.s.len() - 2].clone()) % mi.clone();
        debug!("c{i} = {ci:?}");
        c.push(ci);
    }

    let mut res = c
        .iter()
        .zip(ms)
        .map(|(c_i, m_i)| c_i.clone() * (m.clone() / m_i.clone()))
        .reduce(|a, b| a + b)
        .unwrap()
        % m.clone();
    debug!("The result was {res:?}, but now find the smallest positive integer mod {m:?}");
    while res < R::zero() {
        res = (res + m.clone()) % m.clone();
    }
    res
}
