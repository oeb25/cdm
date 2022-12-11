use crate::euclidean_domain::EuclideanDomain;

#[derive(Debug)]
pub struct ExtendedEuclideanAlgorithm<D> {
    pub r: Vec<D>,
    pub s: Vec<D>,
    pub t: Vec<D>,
    pub q: Vec<D>,
}

impl<D: std::fmt::Debug> std::fmt::Display for ExtendedEuclideanAlgorithm<D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use comfy_table::{
            modifiers::UTF8_ROUND_CORNERS, presets::UTF8_FULL, Attribute, Cell, ContentArrangement,
            Table,
        };

        let mut table = Table::new();
        table
            .load_preset(UTF8_FULL)
            .apply_modifier(UTF8_ROUND_CORNERS)
            .set_content_arrangement(ContentArrangement::Dynamic);

        table.set_header(
            ["i", "q", "r", "s", "t"].map(|t| Cell::new(t).add_attribute(Attribute::Bold)),
        );

        for i in 0..self.r.len() {
            if i == self.r.len() - 1 {
                table.add_row([
                    format!("{i}"),
                    "".to_string(),
                    format!("{:?}", self.r[i]),
                    format!("{:?}", self.s[i]),
                    format!("{:?}", self.t[i]),
                ]);
            } else {
                table.add_row([
                    format!("{i}"),
                    format!("{:?}", self.q[i]),
                    format!("{:?}", self.r[i]),
                    format!("{:?}", self.s[i]),
                    format!("{:?}", self.t[i]),
                ]);
            }
        }

        write!(f, "{table}")
    }
}

// Algorithm 3.6 Traditional Extended Euclidean Algorithm.
pub fn extended_euclidean_algorithm<D: EuclideanDomain + PartialEq + std::fmt::Debug>(
    f: &D,
    g: &D,
) -> ExtendedEuclideanAlgorithm<D> {
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

    ExtendedEuclideanAlgorithm { r, s, t, q }
}

impl<D: EuclideanDomain + PartialEq + std::fmt::Debug> ExtendedEuclideanAlgorithm<D> {
    pub fn gcd(&self) -> &D {
        &self.r[self.r.len() - 2]
    }
    pub fn property(f: D, g: D) -> bool {
        let res = extended_euclidean_algorithm(&f, &g);
        let gcd = res.gcd();

        (f % gcd.clone()).is_zero() && (g % gcd.clone()).is_zero()
    }
}
