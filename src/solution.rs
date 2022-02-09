use std::fmt::{Display, Formatter};
use std::ops::{Index, IndexMut};

pub type Pos = (f64, f64, f64);

#[derive(Clone)]
pub struct VecPos(pub Vec<Pos>);

impl VecPos {
    pub fn with_size(size: usize) -> Self {
        Self(vec![(0f64, 0f64, 0f64); size])
    }
}

impl Display for VecPos {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for (node, position) in self.0.iter().enumerate().skip(1) {
            writeln!(f, "{} {:.9} {:.9} {:.9}", node, position.0, position.1, position.2)?;
        }
        Ok(())
    }
}

impl Index<usize> for VecPos {
    type Output = Pos;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl IndexMut<usize> for VecPos {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}