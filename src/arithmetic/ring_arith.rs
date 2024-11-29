//! This file defines and implements Saber ring elements, specifically ℤ[X]/(X^256 + 1) mod n where
//! n can be any power of two at most 2^16

use crate::{
    consts::RING_DEG,
    ser::{deserialize, serialize},
};

use core::ops::{Add, Mul, Sub};

/// An element of the ring (Z/2^13 Z)[X] / (X^256 + 1)
// The coefficients are in order of ascending powers, i.e., `self.0[0]` is the constant term
#[derive(Debug, Clone, Copy)]
// Auto-derive strict equality for tests. Equality mod q is implemented in ring_eqv in tests below
#[cfg_attr(test, derive(Eq, PartialEq))]
pub struct RingElem(pub(crate) [u16; RING_DEG]);

impl Default for RingElem {
    fn default() -> Self {
        RingElem([0u16; RING_DEG])
    }
}

impl RingElem {
    /// Creates a random ring element
    #[cfg(test)]
    pub(crate) fn rand(rng: &mut impl rand_core::CryptoRngCore) -> Self {
        let modulus = 1 << crate::consts::MODULUS_Q_BITS as u32;

        let mut result = [0; RING_DEG];
        result.iter_mut().for_each(|coeff| {
            let w = rng.next_u32() % modulus;
            *coeff = w as u16;
        });

        RingElem(result)
    }

    /// Deserializes a ring element, treating each coefficient as having only `bits_per_elem` bits.
    /// In Saber terms, this runs BS2POLYk where k = bits_per_elem
    pub(crate) fn from_bytes(bytes: &[u8], bits_per_elem: usize) -> Self {
        assert_eq!(bytes.len(), bits_per_elem * RING_DEG / 8);
        let arr = deserialize(bytes, bits_per_elem);
        RingElem(arr)
    }

    /// Serializes this ring element, treating each coefficient as having only `bits_per_elem`
    /// bits. In Saber terms, this runs POLYk2BS where k = bits_per_elem
    pub(crate) fn to_bytes(self, out_buf: &mut [u8], bits_per_elem: usize) {
        assert_eq!(out_buf.len(), bits_per_elem * RING_DEG / 8);
        serialize(&self.0, out_buf, bits_per_elem)
    }

    // Algorithm 8, ShiftRight
    /// Right-shifts each coefficient by the specified amount, essentially dividing each coeff by a
    /// power of two with rounding
    pub(crate) fn shift_right(&mut self, shift: usize) {
        for coeff in self.0.iter_mut() {
            *coeff >>= shift;
        }
    }

    // Algorithm 7, ShiftLeft
    /// Left-shifts each coefficient by the specified amount, essentially multiplying each coeff by
    /// a power of two, mod 2^16
    pub(crate) fn shift_left(&mut self, shift: usize) {
        for coeff in self.0.iter_mut() {
            *coeff <<= shift;
        }
    }

    /// Adds a given value to all coefficients
    pub(crate) fn wrapping_add_to_all(&mut self, val: u16) {
        for coeff in self.0.iter_mut() {
            *coeff = coeff.wrapping_add(val);
        }
    }
}

impl<'a> Mul for &'a RingElem {
    type Output = RingElem;

    // School book multiplication
    fn mul(self, other: &'a RingElem) -> Self::Output {
        //karatsuba_mul_helper(&self.0, &other.0)
        toom_cook_4way(&self.0, &other.0)
        // Replace the above line with the below line to remove Karatsuba optimization
        //schoolbook_mul_helper(&self.0, &other.0)
    }
}

fn arr_add<const N: usize>(x: &[u16; N], y: &[u16; N]) -> [u16; N] {
    let mut ret = [0u16; N];
    for i in 0..N {
        ret[i] = x[i].wrapping_add(y[i]);
    }
    ret
}

fn arr_sub<const N: usize>(x: &[u16; N], y: &[u16; N]) -> [u16; N] {
    let mut ret = [0u16; N];
    for i in 0..N {
        ret[i] = x[i].wrapping_sub(y[i]);
    }
    ret
}

/// Multiplies a deg-62 polynomial by X^S, where S is at most 64
fn mul63_by_xpow<const S: usize>(p: &[u16; 63]) -> [u16; 127] {
    let mut ret = [0u16; 127];
    for i in 0..63 {
        ret[i + S] = ret[i + S].wrapping_add(p[i]);
    }
    ret
}

/// Schoolbook multiplication of degree-31 polynomials in little-endian order
fn schoolbook32(p: &[u16; 32], q: &[u16; 32]) -> [u16; 63] {
    let mut ret = [0u16; 63];
    // p and q are low-deg enough that there is no wrapping around the ring degree (256)
    for (i, p_coeff) in p.iter().enumerate() {
        for (j, q_coeff) in q.iter().enumerate() {
            let prod = p_coeff.wrapping_mul(*q_coeff);
            ret[i + j] = ret[i + j].wrapping_add(prod);
        }
    }
    return ret;
}

/*
/// Schoolbook multiplication of degree-63 polynomials in little-endian order
fn schoolbook64(p: &[u16; 64], q: &[u16; 64]) -> [u16; 127] {
    let mut result = [0u16; 127];
    // Do all the multiplications
    for (i, p_coeff) in p.iter().enumerate() {
        for (j, q_coeff) in q.iter().enumerate() {
            result[i + j] = result[i + j].wrapping_add(p_coeff.wrapping_mul(*q_coeff));
        }
    }

    result
}
*/

/// Karatsuba multiplication of degree-63 polynomials in little-endian order
fn karatsuba64(p: &[u16; 64], q: &[u16; 64]) -> [u16; 127] {
    const MID: usize = 32;
    const N: usize = 64;

    // If you want to follow along, I used page 4 of these lecture notes
    //     https://cs.dartmouth.edu/~deepc/LecNotes/cs31/lec6.pdf
    // and the Wikipedia page on Karatsuba multiplication
    //     https://en.wikipedia.org/wiki/Karatsuba_algorithm
    // the z_i notation comes from the Wikipedia page

    // Split the inputs into low and high halves
    let pl: &[u16; MID] = &p[..MID].try_into().unwrap();
    let ph: &[u16; MID] = &p[MID..].try_into().unwrap();
    let ql: &[u16; MID] = &q[..MID].try_into().unwrap();
    let qh: &[u16; MID] = &q[MID..].try_into().unwrap();

    // Compute our intermediate products
    let z0 = schoolbook32(pl, ql);
    let z2 = schoolbook32(ph, qh);
    let z3 = schoolbook32(&arr_add(pl, ph), &arr_add(ql, qh));
    let z1 = arr_sub(&arr_sub(&z3, &z2), &z0);

    // Compute z0 + z1*X^mid + z2*X^(2mid)
    let z1 = mul63_by_xpow::<MID>(&z1);
    let z2 = mul63_by_xpow::<N>(&z2);

    let mut padded_z0 = [0u16; 127];
    padded_z0[0..63].copy_from_slice(&z0);
    arr_add(&arr_add(&padded_z0, &z1), &z2)
}

impl<'a> Add for &'a RingElem {
    type Output = RingElem;

    fn add(self, other: &'a RingElem) -> Self::Output {
        RingElem(arr_add(&self.0, &other.0))
    }
}

impl<'a> Sub for &'a RingElem {
    type Output = RingElem;

    fn sub(self, other: &'a RingElem) -> Self::Output {
        RingElem(arr_sub(&self.0, &other.0))
    }
}

const N_SB: usize = RING_DEG / 4;
const N_SB_RES: usize = 2 * N_SB - 1;

#[allow(non_snake_case)]
fn toom_cook_4way(a1: &[u16], b1: &[u16]) -> RingElem {
    let inv3: u32 = 43691;
    let inv9: u32 = 36409;
    let inv15: u32 = 61167;

    let mut aw1 = [0u16; N_SB];
    let mut aw2 = [0u16; N_SB];
    let mut aw3 = [0u16; N_SB];
    let mut aw4 = [0u16; N_SB];
    let mut aw5 = [0u16; N_SB];
    let mut aw6 = [0u16; N_SB];
    let mut aw7 = [0u16; N_SB];

    let mut bw1 = [0u16; N_SB];
    let mut bw2 = [0u16; N_SB];
    let mut bw3 = [0u16; N_SB];
    let mut bw4 = [0u16; N_SB];
    let mut bw5 = [0u16; N_SB];
    let mut bw6 = [0u16; N_SB];
    let mut bw7 = [0u16; N_SB];

    let mut r0;
    let mut r1;
    let mut r2;
    let mut r3;
    let mut r4;
    let mut r5;
    let mut r6;
    let mut r7;

    let (A0, rest) = a1.split_at(N_SB);
    let (A1, rest) = rest.split_at(N_SB);
    let (A2, rest) = rest.split_at(N_SB);
    let (A3, rest) = rest.split_at(N_SB);
    assert!(rest.is_empty());

    let (B0, rest) = b1.split_at(N_SB);
    let (B1, rest) = rest.split_at(N_SB);
    let (B2, rest) = rest.split_at(N_SB);
    let (B3, rest) = rest.split_at(N_SB);
    assert!(rest.is_empty());

    let mut C = [0u16; 2 * RING_DEG];

    // EVALUATION
    for j in 0..N_SB {
        r0 = A0[j];
        r1 = A1[j];
        r2 = A2[j];
        r3 = A3[j];
        r4 = r0.wrapping_add(r2);
        r5 = r1.wrapping_add(r3);
        r6 = r4.wrapping_add(r5);
        r7 = r4.wrapping_sub(r5);
        aw3[j] = r6;
        aw4[j] = r7;
        r4 = (r0 << 2).wrapping_add(r2) << 1;
        r5 = (r1 << 2).wrapping_add(r3);
        r6 = r4.wrapping_add(r5);
        r7 = r4.wrapping_sub(r5);
        aw5[j] = r6;
        aw6[j] = r7;
        r4 = (r3 << 3)
            .wrapping_add(r2 << 2)
            .wrapping_add(r1 << 1)
            .wrapping_add(r0);
        aw2[j] = r4;
        aw7[j] = r0;
        aw1[j] = r3;
    }

    for j in 0..N_SB {
        r0 = B0[j];
        r1 = B1[j];
        r2 = B2[j];
        r3 = B3[j];
        r4 = r0.wrapping_add(r2);
        r5 = r1.wrapping_add(r3);
        r6 = r4.wrapping_add(r5);
        r7 = r4.wrapping_sub(r5);
        bw3[j] = r6;
        bw4[j] = r7;
        r4 = (r0 << 2).wrapping_add(r2) << 1;
        r5 = (r1 << 2).wrapping_add(r3);
        r6 = r4.wrapping_add(r5);
        r7 = r4.wrapping_sub(r5);
        bw5[j] = r6;
        bw6[j] = r7;
        r4 = (r3 << 3)
            .wrapping_add(r2 << 2)
            .wrapping_add(r1 << 1)
            .wrapping_add(r0);
        bw2[j] = r4;
        bw7[j] = r0;
        bw1[j] = r3;
    }

    // MULTIPLICATION

    /*
    let w1 = schoolbook64(&aw1, &bw1);
    let w2 = schoolbook64(&aw2, &bw2);
    let w3 = schoolbook64(&aw3, &bw3);
    let w4 = schoolbook64(&aw4, &bw4);
    let w5 = schoolbook64(&aw5, &bw5);
    let w6 = schoolbook64(&aw6, &bw6);
    let w7 = schoolbook64(&aw7, &bw7);
    */

    let w1 = karatsuba64(&aw1, &bw1);
    let w2 = karatsuba64(&aw2, &bw2);
    let w3 = karatsuba64(&aw3, &bw3);
    let w4 = karatsuba64(&aw4, &bw4);
    let w5 = karatsuba64(&aw5, &bw5);
    let w6 = karatsuba64(&aw6, &bw6);
    let w7 = karatsuba64(&aw7, &bw7);

    // INTERPOLATION
    for i in 0..N_SB_RES {
        r0 = w1[i];
        r1 = w2[i];
        r2 = w3[i];
        r3 = w4[i];
        r4 = w5[i];
        r5 = w6[i];
        r6 = w7[i];

        r1 = r1.wrapping_add(r4);
        r5 = r5.wrapping_sub(r4);
        r3 = r3.wrapping_sub(r2) >> 1;
        r4 = r4.wrapping_sub(r0);
        r4 = r4.wrapping_sub(r6 << 6);
        r4 = (r4 << 1).wrapping_add(r5);
        r2 = r2.wrapping_add(r3);
        r1 = r1.wrapping_sub(r2 << 6).wrapping_sub(r2);
        r2 = r2.wrapping_sub(r6);
        r2 = r2.wrapping_sub(r0);
        r1 = r1.wrapping_add(45u16.wrapping_mul(r2));
        r4 = ((r4.wrapping_sub(r2 << 3) as u32).wrapping_mul(inv3) >> 3) as u16;
        r5 = r5.wrapping_add(r1);
        r1 = ((r1.wrapping_add(r3 << 4) as u32).wrapping_mul(inv9) >> 1) as u16;
        r3 = 0u16.wrapping_sub(r3.wrapping_add(r1));
        r5 = ((30u16.wrapping_mul(r1).wrapping_sub(r5) as u32).wrapping_mul(inv15) >> 2) as u16;
        r2 = r2.wrapping_sub(r4);
        r1 = r1.wrapping_sub(r5);

        C[i] = C[i].wrapping_add(r6);
        C[i + 64] = C[i + 64].wrapping_add(r5);
        C[i + 128] = C[i + 128].wrapping_add(r4);
        C[i + 192] = C[i + 192].wrapping_add(r3);
        C[i + 256] = C[i + 256].wrapping_add(r2);
        C[i + 320] = C[i + 320].wrapping_add(r1);
        C[i + 384] = C[i + 384].wrapping_add(r0);
    }

    let mut res = [0u16; RING_DEG];
    for i in RING_DEG..(2 * RING_DEG) {
        res[i - RING_DEG] = res[i - RING_DEG].wrapping_add(C[i - RING_DEG].wrapping_sub(C[i]));
    }

    RingElem(res)
}

#[cfg(test)]
pub(crate) mod test {
    use super::*;
    use crate::consts::{MODULUS_Q_BITS, RING_DEG};

    use rand::{thread_rng, Rng, RngCore};

    /// Checks that two ring elements are equivalent mod q
    pub(crate) fn ring_eqv(a: RingElem, b: RingElem) -> bool {
        a.0.iter()
            .zip(b.0.iter())
            .all(|(aa, &bb)| aa.wrapping_sub(bb) % (1 << MODULUS_Q_BITS) == 0)
    }

    /// Does schoolbook multiplication of two ring elements
    fn schoolbook_mul_helper(p: &[u16], q: &[u16]) -> RingElem {
        let mut result = RingElem::default();
        // Do all the multiplications
        for (i, p_coeff) in p.iter().enumerate() {
            let mut q_iter = q.iter().enumerate();

            // Do multiplications up until i+j == RING_DEG
            for _ in 0..(RING_DEG - i) {
                let (j, q_coeff) = q_iter.next().expect("there are RING_DEG elems in q_iter");
                let idx = i + j;
                let prod = p_coeff.wrapping_mul(*q_coeff);
                result.0[idx] = result.0[idx].wrapping_add(prod);
            }

            // Once we're past the ring degree, wrap around and multiply by -1. This is
            // because the ring is Z[X]/(X^256 + 1), so X^256 = -1
            for (j, q_coeff) in q_iter {
                let idx = i + j - RING_DEG;
                let prod = p[i].wrapping_mul(*q_coeff);
                result.0[idx] = result.0[idx].wrapping_sub(prod);
            }
        }

        result
    }

    // Checks that a * b == b * a and a + b == b + a for ring elements a, b
    #[test]
    fn commutativity() {
        let mut rng = thread_rng();

        for _ in 0..100 {
            let a = RingElem::rand(&mut rng);
            let b = RingElem::rand(&mut rng);

            let prod_1 = &a * &b;
            let prod_2 = &b * &a;

            let sum_1 = &a + &b;
            let sum_2 = &b + &a;

            assert_eq!(prod_1, prod_2);
            assert_eq!(sum_1, sum_2);
        }
    }

    // Chekcs that ring multiplication is the same as schoolbook multiplciation
    #[test]
    fn mul_correctness() {
        let mut rng = thread_rng();

        for _ in 0..100 {
            let a = RingElem::rand(&mut rng);
            let b = RingElem::rand(&mut rng);

            assert!(ring_eqv(&a * &b, schoolbook_mul_helper(&a.0, &b.0)),);
        }
    }

    // Checks that ring multiplication distributes over addition
    #[test]
    fn distributivity() {
        let mut rng = thread_rng();

        for _ in 0..100 {
            let a = RingElem::rand(&mut rng);
            let b = RingElem::rand(&mut rng);
            let c = RingElem::rand(&mut rng);

            assert!(ring_eqv(&a * &(&b + &c), &(&a * &b) + &(&a * &c)));
        }
    }

    // Tests serialization and deserialization of ring elements
    #[test]
    fn from_bytes() {
        let mut rng = thread_rng();

        // The largest buffer we'll need for the following tests. We make 2 because we need to
        // compare values in some places
        let mut backing_buf1 = [0u8; 16 * RING_DEG / 8];
        let mut backing_buf2 = [0u8; 16 * RING_DEG / 8];

        for _ in 0..1000 {
            // Check that from_bytes matches the reference impl from_bytes for N=2^13,2^10,2^1
            let bits_per_elem = 13;
            let bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
            rng.fill_bytes(bytes);
            assert_eq!(
                reference_impl_from_bytes_mod8192(&bytes),
                RingElem::from_bytes(&bytes, 13)
            );

            // Now check it matches the reference to_bytes impl
            let elem = RingElem::rand(&mut rng);
            let my_bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
            let ref_bytes = &mut backing_buf2[..bits_per_elem * RING_DEG / 8];
            elem.to_bytes(my_bytes, bits_per_elem);
            reference_impl_to_bytes_mod8192(&elem, ref_bytes);
            assert_eq!(my_bytes, ref_bytes);

            let bits_per_elem = 10;
            let bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
            rng.fill_bytes(bytes);
            assert_eq!(
                reference_impl_from_bytes_mod1024(&bytes).0,
                RingElem::from_bytes(&bytes, 10).0,
            );

            let bits_per_elem = 1;
            let bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
            rng.fill_bytes(bytes);
            assert_eq!(
                reference_impl_from_bytes_mod2(&bytes).0,
                RingElem::from_bytes(&bytes, 1).0,
            );

            // Now check it matches the reference to_bytes impl
            let elem = RingElem::rand(&mut rng);
            // The reference impl actually requires that the buffer is zeroed before use
            backing_buf2.fill(0);
            let my_bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
            let ref_bytes = &mut backing_buf2[..bits_per_elem * RING_DEG / 8];
            elem.to_bytes(my_bytes, bits_per_elem);
            reference_impl_to_bytes_mod2(&elem, ref_bytes);
            assert_eq!(my_bytes, ref_bytes);

            // Now check that to_bytes and from_bytes are inverses

            // Pick a random bits_per_elem
            for _ in 0..10 {
                let bits_per_elem = rng.gen_range(1..=13);
                let bitmask = (1 << bits_per_elem) - 1;

                // Generate a random element and make sure none of the values exceed 2^bits_per_elem
                let mut p = RingElem::rand(&mut rng);
                p.0.iter_mut().for_each(|e| *e &= bitmask);

                // Check that a round trip preserves the polynomial
                let p_bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
                p.to_bytes(p_bytes, bits_per_elem);
                assert_eq!(p, RingElem::from_bytes(&p_bytes, bits_per_elem));

                // Now other way around
                let p_bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
                rng.fill_bytes(p_bytes);
                let p = RingElem::from_bytes(&p_bytes, bits_per_elem);
                let new_p_bytes = &mut backing_buf2[..bits_per_elem * RING_DEG / 8];
                p.to_bytes(new_p_bytes, bits_per_elem);
                assert_eq!(p_bytes, new_p_bytes);
            }
        }
    }

    /// A nearly verbatim copy of the C reference impl of BS2POL_N where N = 2^13
    /// https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/pack_unpack.c#L101
    fn reference_impl_from_bytes_mod8192(b: &[u8]) -> RingElem {
        let mut offset_byte;
        let mut offset_data;
        let mut poly = RingElem::default();
        let data = &mut poly.0;

        let b_arr: [u8; 13 * RING_DEG / 8] = b.try_into().unwrap();
        let bytes = b_arr.map(|x| x as u16);

        for j in 0..RING_DEG / 8 {
            offset_byte = 13 * j;
            offset_data = 8 * j;
            data[offset_data] =
                (bytes[offset_byte] & (0xff)) | ((bytes[offset_byte + 1] & 0x1f) << 8);
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

        poly
    }

    /// A nearly verbatim copy of the C reference impl of BS2POL_N where N = 2^10
    /// https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/pack_unpack.c#L134
    fn reference_impl_from_bytes_mod1024(b: &[u8]) -> RingElem {
        let mut offset_byte;
        let mut offset_data;
        let mut poly = RingElem::default();
        let data = &mut poly.0;

        let b_arr: [u8; 10 * RING_DEG / 8] = b.try_into().unwrap();
        let bytes = b_arr.map(|x| x as u16);

        for j in 0..RING_DEG / 4 {
            offset_byte = 5 * j;
            offset_data = 4 * j;
            data[offset_data] =
                (bytes[offset_byte] & (0xff)) | ((bytes[offset_byte + 1] & 0x03) << 8);
            data[offset_data + 1] =
                ((bytes[offset_byte + 1] >> 2) & (0x3f)) | ((bytes[offset_byte + 2] & 0x0f) << 6);
            data[offset_data + 2] =
                ((bytes[offset_byte + 2] >> 4) & (0x0f)) | ((bytes[offset_byte + 3] & 0x3f) << 4);
            data[offset_data + 3] =
                ((bytes[offset_byte + 3] >> 6) & (0x03)) | ((bytes[offset_byte + 4] & 0xff) << 2);
        }

        poly
    }

    /// A nearly verbatim copy of the C reference impl of POL2BS_N where N = 2^13
    /// https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/pack_unpack.c#L78
    fn reference_impl_to_bytes_mod8192(polyn: &RingElem, bytes: &mut [u8]) {
        let mut offset_byte: usize;
        let mut offset_data: usize;
        let data = polyn.0;

        for j in 0..RING_DEG / 8 {
            offset_byte = 13 * j;
            offset_data = 8 * j;
            bytes[offset_byte] = (data[offset_data] & (0xff)) as u8;
            bytes[offset_byte + 1] = ((data[offset_data] >> 8) & 0x1f) as u8
                | ((data[offset_data + 1] & 0x07) << 5) as u8;
            bytes[offset_byte + 2] = ((data[offset_data + 1] >> 3) & 0xff) as u8;
            bytes[offset_byte + 3] = ((data[offset_data + 1] >> 11) & 0x03) as u8
                | ((data[offset_data + 2] & 0x3f) << 2) as u8;
            bytes[offset_byte + 4] = ((data[offset_data + 2] >> 6) & 0x7f) as u8
                | ((data[offset_data + 3] & 0x01) << 7) as u8;
            bytes[offset_byte + 5] = ((data[offset_data + 3] >> 1) & 0xff) as u8;
            bytes[offset_byte + 6] = ((data[offset_data + 3] >> 9) & 0x0f) as u8
                | ((data[offset_data + 4] & 0x0f) << 4) as u8;
            bytes[offset_byte + 7] = ((data[offset_data + 4] >> 4) & 0xff) as u8;
            bytes[offset_byte + 8] = ((data[offset_data + 4] >> 12) & 0x01) as u8
                | ((data[offset_data + 5] & 0x7f) << 1) as u8;
            bytes[offset_byte + 9] = ((data[offset_data + 5] >> 7) & 0x3f) as u8
                | ((data[offset_data + 6] & 0x03) << 6) as u8;
            bytes[offset_byte + 10] = ((data[offset_data + 6] >> 2) & 0xff) as u8;
            bytes[offset_byte + 11] = ((data[offset_data + 6] >> 10) & 0x07) as u8
                | ((data[offset_data + 7] & 0x1f) << 3) as u8;
            bytes[offset_byte + 12] = ((data[offset_data + 7] >> 5) & 0xff) as u8;
        }
    }

    /// A nearly verbatim copy of the C reference impl of BS2POL_N where N = 2
    /// https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/pack_unpack.c#L184
    fn reference_impl_from_bytes_mod2(b: &[u8]) -> RingElem {
        let mut poly = RingElem::default();
        let data = &mut poly.0;

        let b_arr: [u8; 1 * RING_DEG / 8] = b.try_into().unwrap();
        let bytes = b_arr.map(|x| x as u16);

        for j in 0..32 {
            {
                for i in 0..8 {
                    data[j * 8 + i] = (bytes[j] >> i) & 0x01;
                }
            }
        }

        poly
    }

    /// A nearly verbatim copy of the C reference impl of POL2BS_N where N = 2
    /// https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/pack_unpack.c#L196
    fn reference_impl_to_bytes_mod2(polyn: &RingElem, bytes: &mut [u8]) {
        let data = polyn.0;

        for j in 0..32 {
            for i in 0..8 {
                bytes[j] |= ((data[j * 8 + i] & 0x01) << i) as u8;
            }
        }
    }
}
