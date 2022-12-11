use crate::euclidean_domain::{eea, EuclideanDomain};

pub fn chinese_remainder_algorithm<R: EuclideanDomain + PartialOrd>(ms: &[R], v: &[R]) -> R {
    let mut c = vec![];
    let m = ms.iter().cloned().reduce(|m0, m1| m0 * m1).unwrap();
    for i in 0..ms.len() {
        let res = eea(&(m.clone() / ms[i].clone()), &ms[i]);
        println!("{}, {:?}", res, &res.s[res.s.len() - 2]);
        c.push((v[i].clone() * res.s[res.s.len() - 2].clone()) % ms[i].clone());
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
