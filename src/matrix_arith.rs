use crate::ring_arith::{RingElem, RING_DEG};

/// An element of R^{x×y} where R is a [`RingElem`], stored in row-major order
#[derive(Eq, PartialEq, Debug, Clone)]
pub(crate) struct Matrix<const X: usize, const Y: usize>(pub(crate) [[RingElem; X]; Y]);

impl<const X: usize, const Y: usize> Default for Matrix<X, Y> {
    fn default() -> Self {
        Matrix([[RingElem::default(); X]; Y])
    }
}

impl<const X: usize, const Y: usize> Matrix<X, Y> {
    #[cfg(test)]
    pub fn rand(rng: &mut impl rand_core::CryptoRngCore) -> Self {
        let mut mat = Matrix::default();
        for i in 0..Y {
            for j in 0..X {
                mat.0[i][j] = RingElem::rand(rng);
            }
        }
        mat
    }

    /// Applies [`RingElem::shift_right`] to each element in the matrix
    pub(crate) fn shift_right(&mut self, shift: usize) {
        for row in self.0.iter_mut() {
            for elem in row.iter_mut() {
                elem.shift_right(shift)
            }
        }
    }

    /// Multiplies this matrix by the given vector
    pub(crate) fn mul<const Z: usize>(&self, other: &Matrix<Y, Z>) -> Matrix<X, Z> {
        let mut result = Matrix::default();
        for j in 0..X {
            for i in 0..Y {
                for k in 0..Z {
                    let prod = &self.0[i][j] * &other.0[k][i];
                    result.0[k][j] = &result.0[k][j] + &prod;
                }
            }
        }

        result
    }

    /// Multiplies the transpose of this matrix by the given vector
    pub(crate) fn mul_transpose<const Z: usize>(&self, other: &Matrix<X, Z>) -> Matrix<Y, Z> {
        let mut result = Matrix::default();
        for i in 0..Y {
            for j in 0..X {
                for k in 0..Z {
                    let prod = &self.0[i][j] * &other.0[k][j];
                    result.0[k][i] = &result.0[k][i] + &prod;
                }
            }
        }

        result
    }

    /// Adds a given value to all coefficients of all elements of the matrix
    pub(crate) fn wrapping_add_to_all(&mut self, val: u16) {
        for row in self.0.iter_mut() {
            for elem in row.iter_mut() {
                elem.wrapping_add_to_all(val);
            }
        }
    }

    pub(crate) fn to_bytes(&self, out_buf: &mut [u8], bits_per_elem: usize) {
        assert_eq!(out_buf.len(), X * Y * bits_per_elem * RING_DEG / 8);
        let mut chunk_iter = out_buf.chunks_mut(bits_per_elem * RING_DEG / 8);
        for i in 0..Y {
            for j in 0..X {
                let out_chunk = chunk_iter.next().unwrap();
                self.0[i][j].to_bytes(out_chunk, bits_per_elem);
            }
        }
    }

    pub(crate) fn from_bytes(bytes: &[u8], bits_per_elem: usize) -> Self {
        assert_eq!(bytes.len(), X * Y * bits_per_elem * RING_DEG / 8);
        let mut result = Matrix::default();

        let mut chunk_iter = bytes.chunks(bits_per_elem * RING_DEG / 8);
        for i in 0..Y {
            for j in 0..X {
                let chunk = chunk_iter.next().unwrap();
                result.0[i][j] = RingElem::from_bytes(chunk, bits_per_elem);
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

    fn transpose_mat<const X: usize, const Y: usize>(mat: &Matrix<X, Y>) -> Matrix<Y, X> {
        let mut ret = Matrix::default();
        for i in 0..Y {
            for j in 0..X {
                ret.0[j][i] = mat.0[i][j];
            }
        }
        ret
    }

    // Checks that mul and mul_transpose distribute over addition on the RHS
    #[test]
    fn distributivity() {
        const X: usize = 4;
        const Y: usize = 7;

        let mut rng = rand::thread_rng();

        // Test mul_transpose
        let mat = Matrix::<X, Y>::rand(&mut rng);
        let vec1 = Matrix::<X, 1>::rand(&mut rng);
        let vec2 = Matrix::<X, 1>::rand(&mut rng);
        let prod1 = {
            let vec_sum = &vec1 + &vec2;
            mat.mul_transpose(&vec_sum)
        };
        let prod2 = &mat.mul_transpose(&vec1) + &mat.mul_transpose(&vec2);
        assert_eq!(prod1, prod2);

        // Now do the same with mul
        let vec1 = Matrix::<Y, 1>::rand(&mut rng);
        let vec2 = Matrix::<Y, 1>::rand(&mut rng);
        let prod1 = {
            let vec_sum = &vec1 + &vec2;
            mat.mul(&vec_sum)
        };
        let prod2 = &mat.mul(&vec1) + &mat.mul(&vec2);
        assert_eq!(prod1, prod2);
    }

    #[test]
    fn transpose() {
        const X: usize = 4;
        const Y: usize = 7;
        const Z: usize = 13;

        let mut rng = rand::thread_rng();

        // Check that A^T B == (B^T A)^T
        let mat1 = Matrix::<X, Y>::rand(&mut rng);
        let mat2 = Matrix::<X, Z>::rand(&mut rng);
        let prod1 = mat1.mul_transpose(&mat2);
        let prod2 = transpose_mat(&mat2.mul_transpose(&mat1));
        assert_eq!(prod1, prod2);

        // Check that (AB)^T == A^T B^T
        let mat1 = Matrix::<X, Y>::rand(&mut rng);
        let mat2 = Matrix::<Y, Z>::rand(&mut rng);
        let prod1 = transpose_mat(&mat1.mul(&mat2));
        let prod2 = transpose_mat(&mat2).mul(&transpose_mat(&mat1));
        assert_eq!(prod1, prod2);
    }
}
