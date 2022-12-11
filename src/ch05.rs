//! Modular algorithms and interpolation

use crate::euclidean_domain::{eea, EuclideanDomain};

/// ALGORITHM 5.4 Chinese Remainder Algorithm (CRA).
pub fn chinese_remainder_algorithm<R: EuclideanDomain + PartialOrd>(ms: &[R], v: &[R]) -> R {
    let mut c = vec![];
    let m = ms.iter().cloned().reduce(|m0, m1| m0 * m1).unwrap();
    for (mi, vi) in ms.iter().zip(v) {
        let res = eea(&(m.clone() / mi.clone()), &mi);
        println!("{}, {:?}", res, &res.s[res.s.len() - 2]);
        c.push((vi.clone() * res.s[res.s.len() - 2].clone()) % mi.clone());
    }

    let mut res = c
        .iter()
        .zip(ms)
        .map(|(c_i, m_i)| c_i.clone() * (m.clone() / m_i.clone()))
        .reduce(|a, b| a + b)
        .unwrap()
        % m.clone();
    while res < R::zero() {
        res = (res + m.clone()) % m.clone();
    }
    res
}
