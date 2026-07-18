//! This file defines and implements Saber ring elements, specifically ℤ[X]/(X^256 + 1) mod n where
//! n can be any power of two at most 2^16

use crate::{
    consts::RING_DEG,
    ser::{deserialize_generic, serialize},
};

use core::ops::{Add, Mul, Sub};

/// An element of the ring (Z/2^13 Z)[X] / (X^256 + 1)
// The coefficients are in order of ascending powers, i.e., `self.0[0]` is the constant term
#[derive(Eq, PartialEq, Debug, Clone, Copy)]
pub struct RingElem(pub(crate) [u16; RING_DEG]);

impl Default for RingElem {
    fn default() -> Self {
        RingElem([0u16; RING_DEG])
    }
}

impl RingElem {
    /// Creates a random ring element
    #[cfg(test)]
    pub(crate) fn rand(rng: &mut impl rand_core::CryptoRng) -> Self {
        let modulus = 1 << crate::consts::MODULUS_Q_BITS as u32;

        let mut result = [0; RING_DEG];
        result.iter_mut().for_each(|coeff| {
            let w = rng.next_u32() % modulus;
            *coeff = w as u16;
        });

        RingElem(result)
    }

    /// Deserializes a ring element, treating each coefficient as having only `bits_per_elem` bits.
    pub(crate) fn deserialize(bytes: &[u8], bits_per_elem: usize) -> Self {
        assert_eq!(bytes.len(), bits_per_elem * RING_DEG / 8);

        // Specialize based on bits_per_elem. unwraps are okay because of the check aboev
        if bits_per_elem == crate::consts::MODULUS_Q_BITS {
            let arr: &[u8; 13 * RING_DEG / 8] = bytes.try_into().unwrap();
            RingElem(crate::ser::deserialize_13(arr))
        } else if bits_per_elem == crate::consts::MODULUS_P_BITS {
            let arr: &[u8; 10 * RING_DEG / 8] = bytes.try_into().unwrap();
            RingElem(crate::ser::deserialize_10(arr))
        } else {
            RingElem(deserialize_generic(bytes, bits_per_elem))
        }
    }

    /// Serializes this ring element, treating each coefficient as having only `bits_per_elem`
    /// bits. In Saber terms, this runs POLYk2BS where k = bits_per_elem
    pub(crate) fn serialize(&self, out_buf: &mut [u8], bits_per_elem: usize) {
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

    fn mul(self, other: &'a RingElem) -> Self::Output {
        let mut ret = RingElem::default();
        ring_mul_acc(&mut ret, self, other);
        ret
    }
}

// Half the ring degree. We split 256-coefficient polys into two 128-coefficient halves
// for a single level of Karatsuba. Two levels (64×64 base) was benchmarked ~15% slower
// due to increased overhead and less efficient vectorization of shorter inner loops.
const HALF: usize = RING_DEG / 2; // 128

/// Schoolbook multiplication of two 128-coefficient polynomials.
/// The product of two degree-127 polys has degree at most 254, so all 256 output slots suffice.
/// Writes result into `out[0..254]`; `out` must be zeroed on entry.
///
/// Takes fixed-size array references so the compiler knows the exact bounds and can
/// eliminate all bounds checks and vectorize the inner loop.
#[inline(never)] // Benchmarked: keeping this separate lets the compiler vectorize the inner loop better
fn schoolbook_128(out: &mut [u16; RING_DEG], a: &[u16; HALF], b: &[u16; HALF]) {
    // Standard O(n²) schoolbook. The inner loop over b is contiguous in memory, which is
    // cache-friendly. The compiler can hoist a[i] as a loop-invariant broadcast.
    for i in 0..HALF {
        let ai = a[i];
        for j in 0..HALF {
            out[i + j] = out[i + j].wrapping_add(ai.wrapping_mul(b[j]));
        }
    }
}

/// Multiplies two ring elements using one level of Karatsuba, and **accumulates** the product
/// into `acc`. This is the core hot function for Saber's matrix-vector multiplies.
#[allow(clippy::unwrap_used)]
// We use unwrap to split the slices. Once it's stable we should use split_array_ref()
//   https://doc.rust-lang.org/std/primitive.array.html#method.split_array_ref
pub(crate) fn ring_mul_acc(acc: &mut RingElem, a: &RingElem, b: &RingElem) {
    // Convert slices to fixed-size array references for the schoolbook function.
    // These are infallible since we split a RING_DEG array exactly in half.
    let a_lo: &[u16; HALF] = a.0[..HALF].try_into().unwrap();
    let a_hi: &[u16; HALF] = a.0[HALF..].try_into().unwrap();
    let b_lo: &[u16; HALF] = b.0[..HALF].try_into().unwrap();
    let b_hi: &[u16; HALF] = b.0[HALF..].try_into().unwrap();

    // We split each input into low and high 128-coefficient halves:
    //     a = a_lo + a_hi * X^128,   b = b_lo + b_hi * X^128
    // Then use Karatsuba's identity:
    //     a*b = z0 + z1*X^128 + z2*X^256
    // where z0 = a_lo*b_lo, z2 = a_hi*b_hi, z3 = (a_lo+a_hi)*(b_lo+b_hi), z1 = z3 - z0 - z2.
    //
    // Since we work in Z[X]/(X^256 + 1), X^256 = -1, so:
    //     a*b mod (X^256+1) = (z0 - z2) + z1*X^128  mod (X^256+1)
    // And z1*X^128 wraps: coefficients 0..127 of z1 go to positions 128..255,
    // while coefficients 128..255 of z1 wrap to positions 0..127 with a sign flip.

    // Compute the three schoolbook products into flat arrays.
    // Each is a product of two degree-127 polynomials, fitting in 256 coefficients.
    let mut z0 = [0u16; RING_DEG];
    let mut z2 = [0u16; RING_DEG];
    schoolbook_128(&mut z0, a_lo, b_lo);
    schoolbook_128(&mut z2, a_hi, b_hi);

    // Compute (a_lo + a_hi) and (b_lo + b_hi) for the cross term
    let mut a_sum = [0u16; HALF];
    let mut b_sum = [0u16; HALF];
    for i in 0..HALF {
        a_sum[i] = a_lo[i].wrapping_add(a_hi[i]);
        b_sum[i] = b_lo[i].wrapping_add(b_hi[i]);
    }
    let mut z3 = [0u16; RING_DEG];
    schoolbook_128(&mut z3, &a_sum, &b_sum);

    // Accumulate: acc += z0 - z2 + z1*X^128  mod (X^256+1)
    //   where z1 = (z3 - z0 - z2)
    for j in 0..HALF {
        // For acc[j] where j in 0..HALF:
        //   * z0[j] - z2[j] from the direct terms
        //   * -(z3[j+HALF] - z0[j+HALF] - z2[j+HALF]) from z1[j+HALF] wrapping with negation
        let z1_wrap = z3[j + HALF]
            .wrapping_sub(z0[j + HALF])
            .wrapping_sub(z2[j + HALF]);
        acc.0[j] = acc.0[j]
            .wrapping_add(z0[j])
            .wrapping_sub(z2[j])
            .wrapping_sub(z1_wrap);

        // For acc[j] where j in HALF..RING_DEG:
        //   * z0[j] - z2[j] from the direct terms
        //   * +(z3[j-HALF] - z0[j-HALF] - z2[j-HALF]) from z1[j-HALF] (no wrap)
        let z1_direct = z3[j].wrapping_sub(z0[j]).wrapping_sub(z2[j]);
        acc.0[j + HALF] = acc.0[j + HALF]
            .wrapping_add(z0[j + HALF])
            .wrapping_sub(z2[j + HALF])
            .wrapping_add(z1_direct);
    }
}

impl<'a> Add for &'a RingElem {
    type Output = RingElem;

    fn add(self, other: &'a RingElem) -> Self::Output {
        let mut ret = RingElem::default();
        for i in 0..RING_DEG {
            ret.0[i] = self.0[i].wrapping_add(other.0[i]);
        }
        ret
    }
}

impl<'a> Sub for &'a RingElem {
    type Output = RingElem;

    fn sub(self, other: &'a RingElem) -> Self::Output {
        let mut ret = RingElem::default();
        for i in 0..RING_DEG {
            ret.0[i] = self.0[i].wrapping_sub(other.0[i]);
        }
        ret
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::consts::RING_DEG;

    use rand::{rng, Rng, RngCore};

    // Checks that a * b == b * a and a + b == b + a for ring elements a, b
    #[test]
    fn commutativity() {
        let mut rng = rng();

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

    /// Naive schoolbook multiplication directly in Z[X]/(X^256+1) for testing.
    /// This is the simplest correct implementation: O(n^2) with explicit ring reduction.
    fn reference_schoolbook_ring_mul(a: &RingElem, b: &RingElem) -> RingElem {
        let mut result = RingElem::default();
        for i in 0..RING_DEG {
            for j in 0..RING_DEG {
                let prod = a.0[i].wrapping_mul(b.0[j]);
                let idx = i + j;
                if idx < RING_DEG {
                    result.0[idx] = result.0[idx].wrapping_add(prod);
                } else {
                    // X^256 = -1 in our ring, so wrap and negate
                    result.0[idx - RING_DEG] = result.0[idx - RING_DEG].wrapping_sub(prod);
                }
            }
        }
        result
    }

    // Tests that our Karatsuba-based ring_mul_acc matches the naive schoolbook ring multiply
    #[test]
    fn karatsuba_vs_schoolbook() {
        let mut rng = rng();

        for _ in 0..100 {
            let a = RingElem::rand(&mut rng);
            let b = RingElem::rand(&mut rng);

            let reference = reference_schoolbook_ring_mul(&a, &b);
            let optimized = &a * &b;

            assert_eq!(reference, optimized);
        }
    }

    // Tests that ring_mul_acc correctly accumulates into a non-zero buffer
    #[test]
    fn mul_acc_accumulates() {
        let mut rng = rng();

        for _ in 0..50 {
            let a = RingElem::rand(&mut rng);
            let b = RingElem::rand(&mut rng);
            let c = RingElem::rand(&mut rng);
            let d = RingElem::rand(&mut rng);

            // Compute a*b + c*d via ring_mul_acc
            let mut acc = RingElem::default();
            ring_mul_acc(&mut acc, &a, &b);
            ring_mul_acc(&mut acc, &c, &d);

            // Compute the same thing via separate multiplies + add
            let expected = &(&a * &b) + &(&c * &d);

            assert_eq!(acc, expected);
        }
    }

    // Tests serialization and deserialization of ring elements
    #[test]
    fn deserialize() {
        let mut rng = rng();

        // The largest buffer we'll need for the following tests. We make 2 because we need to
        // compare values in some places
        let mut backing_buf1 = [0u8; 16 * RING_DEG / 8];
        let mut backing_buf2 = [0u8; 16 * RING_DEG / 8];

        for _ in 0..1000 {
            // Check that deserialize matches the reference impl deserialize for N=2^13,2^10,2^1
            let bits_per_elem = 13;
            let bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
            rng.fill_bytes(bytes);
            assert_eq!(
                saber_ref_from_bytes_mod8192(&bytes),
                RingElem::deserialize(&bytes, 13)
            );

            // Now check it matches the reference to_bytes impl
            let elem = RingElem::rand(&mut rng);
            let my_bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
            let ref_bytes = &mut backing_buf2[..bits_per_elem * RING_DEG / 8];
            elem.serialize(my_bytes, bits_per_elem);
            reference_impl_to_bytes_mod8192(&elem, ref_bytes);
            assert_eq!(my_bytes, ref_bytes);

            let bits_per_elem = 10;
            let bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
            rng.fill_bytes(bytes);
            assert_eq!(
                saber_ref_from_bytes_mod1024(&bytes).0,
                RingElem::deserialize(&bytes, 10).0,
            );

            let bits_per_elem = 1;
            let bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
            rng.fill_bytes(bytes);
            assert_eq!(
                saber_ref_from_bytes_mod2(&bytes).0,
                RingElem::deserialize(&bytes, 1).0,
            );

            // Now check it matches the reference to_bytes impl
            let elem = RingElem::rand(&mut rng);
            // The reference impl actually requires that the buffer is zeroed before use
            backing_buf2.fill(0);
            let my_bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
            let ref_bytes = &mut backing_buf2[..bits_per_elem * RING_DEG / 8];
            elem.serialize(my_bytes, bits_per_elem);
            saber_ref_to_bytes_mod2(&elem, ref_bytes);
            assert_eq!(my_bytes, ref_bytes);

            // Now check that to_bytes and from_bytes are inverses

            // Pick a random bits_per_elem
            for _ in 0..10 {
                let bits_per_elem = rng.random_range(1..=13);
                let bitmask = (1 << bits_per_elem) - 1;

                // Generate a random element and make sure none of the values exceed 2^bits_per_elem
                let mut p = RingElem::rand(&mut rng);
                p.0.iter_mut().for_each(|e| *e &= bitmask);

                // Check that a round trip preserves the polynomial
                let p_bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
                p.serialize(p_bytes, bits_per_elem);
                assert_eq!(p, RingElem::deserialize(&p_bytes, bits_per_elem));

                // Now other way around
                let p_bytes = &mut backing_buf1[..bits_per_elem * RING_DEG / 8];
                rng.fill_bytes(p_bytes);
                let p = RingElem::deserialize(&p_bytes, bits_per_elem);
                let new_p_bytes = &mut backing_buf2[..bits_per_elem * RING_DEG / 8];
                p.serialize(new_p_bytes, bits_per_elem);
                assert_eq!(p_bytes, new_p_bytes);
            }
        }
    }

    /// A nearly verbatim copy of the C reference impl of BS2POL_N where N = 2^13
    /// https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/pack_unpack.c#L101
    fn saber_ref_from_bytes_mod8192(b: &[u8]) -> RingElem {
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
    fn saber_ref_from_bytes_mod1024(b: &[u8]) -> RingElem {
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
    fn saber_ref_from_bytes_mod2(b: &[u8]) -> RingElem {
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
    fn saber_ref_to_bytes_mod2(polyn: &RingElem, bytes: &mut [u8]) {
        let data = polyn.0;

        for j in 0..32 {
            for i in 0..8 {
                bytes[j] |= ((data[j * 8 + i] & 0x01) << i) as u8;
            }
        }
    }
}
