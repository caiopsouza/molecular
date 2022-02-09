use std::fmt::{Debug, Formatter};
use std::ops::{Index, IndexMut};

#[derive(Default, Clone)]
pub struct Torsion([f64; 12]);

impl Torsion {
    // Assumes the rest of the matrix is already zeroed.
    pub(crate) fn eye() -> Self {
        let mut res = Self::default();

        res[(0, 0)] = 1f64;
        res[(1, 1)] = 1f64;
        res[(2, 2)] = 1f64;

        res
    }

    fn product_line(&mut self, one: &Self, other: &Self, line: usize) {
        self[(line, 0)] = one[(line, 0)] * other[(0, 0)] + one[(line, 1)] * other[(1, 0)] + one[(line, 2)] * other[(2, 0)];
        self[(line, 1)] = one[(line, 0)] * other[(0, 1)] + one[(line, 1)] * other[(1, 1)] + one[(line, 2)] * other[(2, 1)];
        self[(line, 2)] = one[(line, 0)] * other[(0, 2)] + one[(line, 1)] * other[(1, 2)] + one[(line, 2)] * other[(2, 2)];
        self[(line, 3)] = one[(line, 0)] * other[(0, 3)] + one[(line, 1)] * other[(1, 3)] + one[(line, 2)] * other[(2, 3)] + one[(line, 3)];
    }

    pub fn product(&self, other: &Self) -> Self {
        let mut res = Self::default();

        res.product_line(self, other, 0);
        res.product_line(self, other, 1);
        res.product_line(self, other, 2);

        res
    }
}

impl Debug for Torsion {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "[[{}, {}, {}, {}]\n [{}, {}, {}, {}]\n [{}, {}, {}, {}]\n [0, 0, 0, 1]]",
               self[(0, 0)],
               self[(0, 1)],
               self[(0, 2)],
               self[(0, 3)],
               self[(1, 0)],
               self[(1, 1)],
               self[(1, 2)],
               self[(1, 3)],
               self[(2, 0)],
               self[(2, 1)],
               self[(2, 2)],
               self[(2, 3)],
        )
    }
}

impl Index<(usize, usize)> for Torsion {
    type Output = f64;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.0[4 * index.0 + index.1]
    }
}

impl IndexMut<(usize, usize)> for Torsion {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.0[4 * index.0 + index.1]
    }
}