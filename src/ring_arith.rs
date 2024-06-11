use core::ops::{Add, Mul};

use rand_core::CryptoRngCore;

/// The modulus of our base ring Z/2^13 Z
const MODULUS: u16 = 1 << 13;

/// The value of -1 in Z/2^13 Z
const NEG_ONE: u16 = MODULUS - 1;

/// The degree of the polynomial ring over Z/2^13 Z
pub(crate) const RING_DEG: usize = 256;

const KARATSUBA_THRESHOLD: usize = 128;

/// An element of the ring (Z/2^13 Z)[X] / (X^256 + 1)
// The coefficients are in order of ascending powers, i.e., `self.0[0]` is the constant term
#[derive(Eq, PartialEq, Debug, Clone, Copy)]
pub struct RingElem([u16; RING_DEG]);

impl Default for RingElem {
    fn default() -> Self {
        RingElem([0u16; RING_DEG])
    }
}

impl RingElem {
    /// Creates a random ring element
    pub fn rand(rng: &mut impl CryptoRngCore) -> Self {
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
        schoolbook_mul_helper(&self.0, &other.0)
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

/// Multiplies the given ring element by X^pow. In our representation, this means shifting the
/// coefficients of p to the right, and multiplying by -1 when they wrap around
fn mul_by_xpow(p: &RingElem, shift: usize) -> RingElem {
    let mut ret = RingElem::default();
    for i in 0..RING_DEG {
        let is_neg = ((i + shift) / RING_DEG) % 2 == 1;
        let idx = (i + shift) % RING_DEG;
        if is_neg {
            ret.0[idx] = ret.0[idx].wrapping_sub(p.0[i]);
        } else {
            ret.0[idx] = ret.0[idx].wrapping_add(p.0[i]);
        }
    }
    ret
}

// Returns p*q and the size (deg+1) of the resulting polyn
fn kara_mul_helper(p: &[u16], q: &[u16]) -> RingElem {
    assert_eq!(p.len(), q.len());
    let n = p.len();

    // Eventually, we have few enough terms that we should just do schoolbook multiplication
    if n == KARATSUBA_THRESHOLD {
        let mut ret = RingElem::default();
        // n is small enough that there is no wrapping around the ring degree (256)
        for i in 0..n {
            for j in 0..n {
                let prod = p[i].wrapping_mul(q[j]);
                ret.0[i + j] = ret.0[i + j].wrapping_add(prod);
            }
        }
        return ret;
    }

    let mid = n / 2;
    let pl = &p[..mid];
    let ph = &p[mid..];
    let ql = &q[..mid];
    let qh = &q[mid..];

    let z0 = kara_mul_helper(pl, ql);
    let z2 = kara_mul_helper(ph, qh);
    let z3 = kara_mul_helper(&poly_add(pl, ph).0[..mid], &poly_add(ql, qh).0[..mid]);
    let z1 = poly_sub(&poly_sub(&z3.0, &z2.0).0, &z0.0);

    // Compute z0 + z1*X^mid + z2*X^(2mid)
    let z1 = mul_by_xpow(&z1, mid);
    let z2 = mul_by_xpow(&z2, 2 * mid);
    let res = poly_add(&poly_add(&z0.0, &z1.0).0, &z2.0);

    res
}

fn schoolbook_mul_helper(p: &[u16], q: &[u16]) -> RingElem {
    let mut result = RingElem::default();
    // Do all the multiplications
    for i in 0..RING_DEG {
        for j in 0..(RING_DEG - i) {
            let idx = i + j;
            let prod = p[i].wrapping_mul(q[j]);
            result.0[idx] = result.0[idx].wrapping_add(prod);
        }
        for j in (RING_DEG - i)..RING_DEG {
            let prod = p[i].wrapping_mul(q[j]);
            let idx = i + j - RING_DEG;
            result.0[idx] = result.0[idx].wrapping_sub(prod);
        }
    }

    result
}

impl RingElem {
    pub fn kara_mul(&self, other: &RingElem) -> RingElem {
        kara_mul_helper(&self.0, &other.0)
    }

    pub fn schoolbook_mul(&self, other: &RingElem) -> RingElem {
        schoolbook_mul_helper(&self.0, &other.0)
    }

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

        for _ in 0..2 {
            let a = RingElem::rand(&mut rng);
            let b = RingElem::rand(&mut rng);

            let mut kara_prod = kara_mul_helper(&a.0, &b.0);
            let mut schoolbook_prod = &a * &b;
            kara_prod.reduce();
            schoolbook_prod.reduce();

            assert_eq!(kara_prod, schoolbook_prod);
        }
    }
}
