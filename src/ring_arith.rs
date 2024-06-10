use core::ops::{Add, Mul};

use rand_core::CryptoRngCore;

/// The modulus of our base ring Z/2^13 Z
const MODULUS: u16 = 1 << 13;

/// The value of -1 in Z/2^13 Z
const NEG_ONE: u16 = MODULUS - 1;

/// The degree of the polynomial ring over Z/2^13 Z
pub(crate) const RING_DEG: usize = 256;

/// An element of the ring (Z/2^13 Z)[X] / (X^256 + 1)
// The coefficients are in order of ascending powers, i.e., `self.0[0]` is the constant term
#[derive(Eq, PartialEq, Debug, Clone, Copy)]
pub(crate) struct RingElem([u16; RING_DEG]);

impl Default for RingElem {
    fn default() -> Self {
        RingElem([0u16; RING_DEG])
    }
}

impl RingElem {
    /// Creates a random ring element
    pub(crate) fn rand(rng: &mut impl CryptoRngCore) -> Self {
        let mut result = [0; RING_DEG];
        for i in 0..RING_DEG {
            let coeff = rng.next_u32() % MODULUS as u32;
            result[i] = coeff as u16;
        }
        RingElem(result)
    }
}

impl<'a> Mul for &'a RingElem {
    type Output = RingElem;

    // School book multiplication
    fn mul(self, other: &'a RingElem) -> Self::Output {
        let mut result = [0u16; RING_DEG];
        // Do all the multiplications
        for i in 0..RING_DEG {
            for j in 0..RING_DEG {
                // We can multiply and add with wrapping because everything is mod a power of 2
                let mut prod = self.0[i].wrapping_mul(other.0[j]);
                let idx = (i + j) % RING_DEG;
                let idx_is_past_deg = (i + j) >= RING_DEG;
                // The coeff a·X^{256 + i} equals -a·X^i since X^256 + 1 == 0 in our polyn ring
                if idx_is_past_deg {
                    prod = prod.wrapping_mul(NEG_ONE);
                }
                result[idx as usize] = result[idx as usize].wrapping_add(prod);
            }
        }

        RingElem(result)
    }
}

fn poly_add(x: &[u16], y: &[u16]) -> RingElem {
    let mut ret = RingElem::default();
    let outlen = core::cmp::max(x.len(), y.len());
    for i in 0..outlen {
        ret.0[i] = x.get(i).unwrap_or(&0).wrapping_add(*y.get(i).unwrap_or(&0))
    }
    ret
}

fn poly_sub(x: &[u16], y: &[u16]) -> RingElem {
    let mut ret = RingElem::default();
    let outlen = core::cmp::max(x.len(), y.len());
    for i in 0..outlen {
        ret.0[i] = x.get(i).unwrap_or(&0).wrapping_sub(*y.get(i).unwrap_or(&0))
    }
    ret
}

fn shift_right(x: &RingElem, shift: usize) -> RingElem {
    let mut ret = RingElem::default();
    for i in 0..RING_DEG {
        let is_neg = ((i + shift) / RING_DEG) % 2 == 1;
        let coeff = if is_neg { NEG_ONE } else { 1u16 };
        ret.0[(i + shift) % RING_DEG] = coeff.wrapping_mul(x.0[i]);
    }
    ret
}

// Returns x*y and the size (deg+1) of the resulting polyn
fn kara_mul(x: &[u16], y: &[u16]) -> RingElem {
    assert_eq!(x.len(), y.len());
    let n = x.len();

    if n == 1 {
        let mut ret = RingElem::default();
        ret.0[0] = x[0].wrapping_mul(y[0]);
        return ret;
    }

    let mid = n / 2;
    let xl = &x[..mid];
    let xh = &x[mid..];
    let yl = &y[..mid];
    let yh = &y[mid..];

    let z0 = kara_mul(xl, yl);
    let z2 = kara_mul(xh, yh);
    let z3 = kara_mul(&poly_add(xl, xh).0[..mid], &poly_add(yl, yh).0[..mid]);
    let z1 = poly_sub(&poly_sub(&z3.0, &z2.0).0, &z0.0);

    // Compute z0 + z1*X^mid + z2*X^(2mid)
    let z1 = shift_right(&z1, mid);
    let z2 = shift_right(&z2, 2 * mid);
    let res = poly_add(&poly_add(&z0.0, &z1.0).0, &z2.0);

    res
}

impl RingElem {
    fn reduce(&mut self) {
        // Now reduce everything
        for i in 0..RING_DEG {
            self.0[i] %= MODULUS;
        }
    }
}

impl<'a> Add for &'a RingElem {
    type Output = RingElem;

    fn add(self, other: &'a RingElem) -> Self::Output {
        let mut result = [0u16; RING_DEG];
        self.0
            .iter()
            .zip(other.0.iter())
            .enumerate()
            .for_each(|(i, (a, b))| result[i] = a.wrapping_add(*b));

        RingElem(result)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use rand::thread_rng;

    // Checks that a * b == b * a and a + b == b + a for ring elements a, b
    #[test]
    fn commutativity() {
        let mut rng = thread_rng();

        for _ in 0..100 {
            let a = RingElem::rand(&mut rng);
            let b = RingElem::rand(&mut rng);

            let mut prod_1 = &a * &b;
            let mut prod_2 = &b * &a;
            prod_1.reduce();
            prod_2.reduce();

            let mut sum_1 = &a + &b;
            let mut sum_2 = &b + &a;
            sum_1.reduce();
            sum_2.reduce();

            assert_eq!(prod_1, prod_2);
            assert_eq!(sum_1, sum_2);
        }
    }

    #[test]
    fn karatsuba() {
        let mut rng = thread_rng();

        for _ in 0..100 {
            let a = RingElem::rand(&mut rng);
            let b = RingElem::rand(&mut rng);

            let k = kara_mul(&a.0, &b.0);
            assert_eq!(k, &a * &b);
        }
    }
}
