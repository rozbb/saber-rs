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

    fn from_bytes_mod8192(bytes: &[u8]) -> Self {
        assert_eq!(bytes.len(), 13 * 256 / 8);

        let mut poly = RingElem::default();
        for (bs, cs) in bytes.chunks_exact(13).zip(poly.0.chunks_exact_mut(8)) {
            cs[0] = (bs[0] as u16) | ((bs[1] as u16) << 8);
            cs[1] = ((bs[1] as u16) >> 5) | ((bs[2] as u16) << 3) | ((bs[3] as u16) << 11);
            cs[2] = ((bs[3] as u16) >> 2) | ((bs[4] as u16) << 6);
            cs[3] = ((bs[4] as u16) >> 7) | ((bs[5] as u16) << 1) | ((bs[6] as u16) << 9);
            cs[4] = ((bs[6] as u16) >> 4) | ((bs[7] as u16) << 4) | ((bs[8] as u16) << 12);
            cs[5] = ((bs[8] as u16) >> 1) | ((bs[9] as u16) << 7);
            cs[6] = ((bs[9] as u16) >> 6) | ((bs[10] as u16) << 2) | ((bs[11] as u16) << 10);
            cs[7] = ((bs[11] as u16) >> 3) | ((bs[12] as u16) << 5);
        }
        poly.reduce();
        poly
    }

    fn from_bytes(bytes: &[u8], bits_per_elem: usize) -> Self {
        assert_eq!(bytes.len(), bits_per_elem * 256 / 8);
        let mut p = RingElem::default();

        // Accumulate all the bits into p
        let mut bit_idx = 0;
        while bit_idx < bits_per_elem * 256 {
            let byte_idx = bit_idx / 8;
            let elem_idx = bit_idx / bits_per_elem;
            let bit_in_byte = bit_idx % 8;
            let bit_in_elem = bit_idx % bits_per_elem;

            p.0[elem_idx] |= ((bytes[byte_idx] as u16) >> bit_in_byte) << bit_in_elem;
            let just_accumulated = core::cmp::min(8 - bit_in_byte, bits_per_elem - bit_in_elem);
            bit_idx += just_accumulated
        }

        p.reduce();
        p
    }

    fn to_bytes(self, out_buf: &mut [u8], bits_per_elem: usize) {
        assert_eq!(out_buf.len(), bits_per_elem * 256 / 8);

        // Write all the bits into the given bytestring
        let mut bit_idx = 0;
        while bit_idx < bits_per_elem * 256 {
            let byte_idx = bit_idx / 8;
            let elem_idx = bit_idx / bits_per_elem;
            let bit_in_byte = bit_idx % 8;
            let bit_in_elem = bit_idx % bits_per_elem;

            out_buf[byte_idx] |= ((self.0[elem_idx] >> bit_in_elem) as u8) << bit_in_byte;
            let just_wrote = core::cmp::min(8 - bit_in_byte, bits_per_elem - bit_in_elem);
            bit_idx += just_wrote
        }
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
fn karatsuba_mul_helper(p: &[u16], q: &[u16]) -> RingElem {
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

    let z0 = karatsuba_mul_helper(pl, ql);
    let z2 = karatsuba_mul_helper(ph, qh);
    let z3 = karatsuba_mul_helper(&poly_add(pl, ph).0[..mid], &poly_add(ql, qh).0[..mid]);
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
    pub fn karatsuba_mul(&self, other: &RingElem) -> RingElem {
        karatsuba_mul_helper(&self.0, &other.0)
    }

    pub fn schoolbook_mul(&self, other: &RingElem) -> RingElem {
        schoolbook_mul_helper(&self.0, &other.0)
    }

    fn reduce(&mut self) {
        // We can do w & MODULUS- 1 instead of w % MODULUS because both simply leave the bottom 13
        // bits
        self.0.iter_mut().for_each(|w| *w &= MODULUS - 1)
    }
}

impl<'a> Add for &'a RingElem {
    type Output = RingElem;

    fn add(self, other: &'a RingElem) -> Self::Output {
        poly_add(&self.0, &other.0)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use rand::{thread_rng, Rng, RngCore};

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

            let mut kara_prod = karatsuba_mul_helper(&a.0, &b.0);
            let mut schoolbook_prod = &a * &b;
            kara_prod.reduce();
            schoolbook_prod.reduce();

            assert_eq!(kara_prod, schoolbook_prod);
        }
    }

    /// A nearly verbatim copy of the C reference impl of BS2POL_N where N = 2^13
    /// https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/pack_unpack.c#L101
    fn refernce_impl_from_bytes_mod8192(b: &[u8]) -> RingElem {
        let mut offset_byte;
        let mut offset_data;
        let mut poly = RingElem::default();
        let data = &mut poly.0;

        let bytes: Vec<u16> = b.iter().map(|&x| x as u16).collect();

        for j in 0..RING_DEG / 8 {
            offset_byte = 13 * j;
            offset_data = 8 * j;
            data[offset_data + 0] =
                (bytes[offset_byte + 0] & (0xff)) | ((bytes[offset_byte + 1] & 0x1f) << 8);
            data[offset_data + 1] = (bytes[offset_byte + 1] >> 5 & (0x07))
                | ((bytes[offset_byte + 2] & 0xff) << 3)
                | ((bytes[offset_byte + 3] & 0x03) << 11);
            data[offset_data + 2] =
                (bytes[offset_byte + 3] >> 2 & (0x3f)) | ((bytes[offset_byte + 4] & 0x7f) << 6);
            data[offset_data + 3] = (bytes[offset_byte + 4] >> 7 & (0x01))
                | ((bytes[offset_byte + 5] & 0xff) << 1)
                | ((bytes[offset_byte + 6] & 0x0f) << 9);
            data[offset_data + 4] = (bytes[offset_byte + 6] >> 4 & (0x0f))
                | ((bytes[offset_byte + 7] & 0xff) << 4)
                | ((bytes[offset_byte + 8] & 0x01) << 12);
            data[offset_data + 5] =
                (bytes[offset_byte + 8] >> 1 & (0x7f)) | ((bytes[offset_byte + 9] & 0x3f) << 7);
            data[offset_data + 6] = (bytes[offset_byte + 9] >> 6 & (0x03))
                | ((bytes[offset_byte + 10] & 0xff) << 2)
                | ((bytes[offset_byte + 11] & 0x07) << 10);
            data[offset_data + 7] =
                (bytes[offset_byte + 11] >> 3 & (0x1f)) | ((bytes[offset_byte + 12] & 0xff) << 5);
        }

        poly.reduce();
        poly
    }

    #[test]
    fn from_bytes() {
        for _ in 0..100 {
            let bits_per_elem = 13;
            let mut rng = thread_rng();
            let mut bytes = vec![0u8; bits_per_elem * 256 / 8];
            rng.fill_bytes(&mut bytes);

            // Check that from_bytes matches the reference impl from_bytes for N=2^13
            assert_eq!(
                refernce_impl_from_bytes_mod8192(&bytes),
                RingElem::from_bytes(&bytes, 13)
            );

            // Now check that to_bytes and from_bytes are inverses. Pick a random bits_per_elem
            let bits_per_elem = rng.gen_range(1..=16);
            let x = RingElem::rand(&mut rng);
            let mut x_bytes = vec![0u8; bits_per_elem * 256 / 8];
            x.to_bytes(&mut x_bytes, dbg!(bits_per_elem));
            assert_eq!(x, RingElem::from_bytes(&x_bytes, bits_per_elem));
        }
    }
}
