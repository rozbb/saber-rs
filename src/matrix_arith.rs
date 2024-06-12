use crate::ring_arith::RingElem;

/// An element of R^{x×y} where R is a [`RingElem`], stored in row-major order
#[derive(Eq, PartialEq, Debug, Clone)]
pub(crate) struct Matrix<const X: usize, const Y: usize>(pub(crate) [[RingElem; X]; Y]);

impl<const X: usize, const Y: usize> Default for Matrix<X, Y> {
    fn default() -> Self {
        Matrix([[RingElem::default(); X]; Y])
    }
}

impl<const X: usize, const Y: usize> Matrix<X, Y> {
    /// Applies [`RingElem::shift_right`] to each element in the matrix
    fn shift_right(&mut self, shift: usize) {
        for mut row in self.0 {
            for mut elem in row {
                elem.shift_right(shift)
            }
        }
    }

    /// Multiplies the transpose of this matrix by the given vector
    fn mul_transpose(&self, other: &Matrix<X, 1>) -> Matrix<Y, 1> {
        let mut result = Matrix::default();
        for i in 0..Y {
            for j in 0..X {
                let prod = &self.0[i][j] * &other.0[0][j];
                result.0[0][i] = &result.0[0][i] + &prod;
            }
        }

        result
    }
}

impl<'a, const X: usize, const Y: usize> core::ops::Add<&'a Matrix<X, Y>> for &'a Matrix<X, Y> {
    type Output = Matrix<X, Y>;

    fn add(self, other: &'a Matrix<X, Y>) -> Self::Output {
        let mut result = Matrix::default();
        for i in 0..Y {
            for j in 0..X {
                result.0[i][j] = &self.0[i][j] + &other.0[i][j];
            }
        }

        result
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use rand::RngCore;

    use crate::gen::{gen_matrix_from_seed, gen_secret_from_seed};

    // Checks that mul_transpose distributes over addition on the RHS
    #[test]
    fn distributivity() {
        const L: usize = 4;
        const MU: usize = 10;

        let mut rng = rand::thread_rng();

        let mut mat_seed = [0u8; 32];
        let mut vec1_seed = [0u8; 32];
        let mut vec2_seed = [0u8; 32];
        rng.fill_bytes(&mut mat_seed);
        rng.fill_bytes(&mut vec1_seed);
        rng.fill_bytes(&mut vec2_seed);

        let mat = gen_matrix_from_seed::<L>(mat_seed);
        let vec1 = gen_secret_from_seed::<L, MU>(vec1_seed);
        let vec2 = gen_secret_from_seed::<L, MU>(vec2_seed);

        let prod1 = {
            let vec_sum = &vec1 + &vec2;
            mat.mul_transpose(&vec_sum)
        };
        let prod2 = &mat.mul_transpose(&vec1) + &mat.mul_transpose(&vec2);

        assert_eq!(prod1, prod2);
    }
}
