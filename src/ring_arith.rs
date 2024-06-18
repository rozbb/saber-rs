use core::ops::{Add, Mul, Sub};

use rand_core::CryptoRngCore;

/// The modulus q of our base ring Z/2^13 Z
pub(crate) const MODULUS_Q: u16 = 1 << 13;

/// The modulus p of some of our public values
pub(crate) const MODULUS_P: u16 = 1 << 10;

/// The bitlength of the modulus q
pub(crate) const MODULUS_Q_BITS: usize = 13;

/// The bitlength of the modulus p
pub const MODULUS_P_BITS: usize = 10;

/// The degree of the polynomial ring over Z/2^13 Z
pub const RING_DEG: usize = 256;

/// The degree (-1) of a polynomial at which point we revert to schoolbook multiplication
const KARATSUBA_THRESHOLD: usize = 128;

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
    pub fn rand(rng: &mut impl CryptoRngCore) -> Self {
        let mut result = [0; RING_DEG];
        for i in 0..RING_DEG {
            let coeff = rng.next_u32() % MODULUS_Q as u32;
            result[i] = coeff as u16;
        }
        RingElem(result)
    }

    /// Deserializes a ring element, treating each element as having only `bits_per_elem` bits.
    /// In Saber terms, this runs BS2POLYk where k = bits_per_elem
    pub(crate) fn from_bytes(bytes: &[u8], bits_per_elem: usize) -> Self {
        assert_eq!(bytes.len(), bits_per_elem * RING_DEG / 8);
        let arr = deserialize(bytes, bits_per_elem);
        RingElem(arr)
    }

    /// Serializes this ring element, treating each element as having only `bits_per_elem` bits.
    /// In Saber terms, this runs POLYk2BS where k = bits_per_elem
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

    /// Reduces the coefficients of this ring element to their canonical representative modulo
    /// 2^modbits
    pub(crate) fn reduce_mod_2pow(&mut self, mod_bits: usize) {
        let bitmask = (1 << mod_bits) - 1;
        self.0.iter_mut().for_each(|w| *w &= bitmask);
    }
}

impl<'a> Mul for &'a RingElem {
    type Output = RingElem;

    // School book multiplication
    fn mul(self, other: &'a RingElem) -> Self::Output {
        schoolbook_mul_helper(&self.0, &other.0)
        //karatsuba_mul_helper(&self.0, &other.0)
    }
}

/// Adds the two inputs, treating them as the lower coefficients of a RingElem
fn poly_add(x: &[u16], y: &[u16]) -> RingElem {
    let mut ret = RingElem::default();
    let outlen = core::cmp::max(x.len(), y.len());
    for i in 0..outlen {
        let lhs = x.get(i).unwrap_or(&0);
        let rhs = y.get(i).unwrap_or(&0);
        ret.0[i] = lhs.wrapping_add(*rhs);
    }
    ret
}

/// Subtracts the two inputs, treating them as the lower coefficients of a RingElem
fn poly_sub(x: &[u16], y: &[u16]) -> RingElem {
    let mut ret = RingElem::default();
    let outlen = core::cmp::max(x.len(), y.len());
    for i in 0..outlen {
        let lhs = x.get(i).unwrap_or(&0);
        let rhs = y.get(i).unwrap_or(&0);
        ret.0[i] = lhs.wrapping_sub(*rhs);
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

/// Returns p*q as ring elements. p and q MUST be the same length
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

    // Split the inputs into low and high halves
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
            let idx = i + j - RING_DEG;
            let prod = p[i].wrapping_mul(q[j]);
            result.0[idx] = result.0[idx].wrapping_sub(prod);
        }
    }

    result
}

/// Deserializes the given bitstring into a u16 array. Every element of the array has
/// `bits_per_elem` bits (must be ≤ 16), encoded in the lower bits of the word.
pub(crate) fn deserialize<const N: usize>(bytes: &[u8], bits_per_elem: usize) -> [u16; N] {
    assert_eq!(bytes.len(), bits_per_elem * N / 8);

    // We only want the lower bits_per_elem bits to be set in any elem of our output
    let bitmask = (1 << bits_per_elem) - 1;

    // Accumulate all the bits into p
    let mut p = [0u16; N];
    let mut bit_idx = 0;
    while bit_idx < bits_per_elem * N {
        let byte_idx = bit_idx / 8;
        let elem_idx = bit_idx / bits_per_elem;
        let bit_in_byte = bit_idx % 8;
        let bit_in_elem = bit_idx % bits_per_elem;

        // Shift the byte we're reading so the first bit we want is the lowest bit
        let data_to_read = bytes[byte_idx] as u16 >> bit_in_byte;
        // OR the byte into our element, shifitng to align with the first unused bit in the elem
        p[elem_idx] |= data_to_read << bit_in_elem;

        // The above might set some high bits we don't want. Clear them now.
        p[elem_idx] &= bitmask;

        // We have read either: 1) however many bits were remaining in the byte we were
        // reading, or 2) however many unused bits remained in the current element we were
        // writing to. Whichever is smaller.
        let just_read = core::cmp::min(8 - bit_in_byte, bits_per_elem - bit_in_elem);
        bit_idx += just_read;
    }

    p
}

/// Serializes the given u16 array into a bitstring. Every element of the array has `bits_per_elem`
/// bits (must be ≤ 16), encoded in the lower bits of the word.
pub(crate) fn serialize(data: &[u16], out_buf: &mut [u8], bits_per_elem: usize) {
    assert_eq!(out_buf.len(), bits_per_elem * data.len() / 8);

    // We only want to write the lower bits_per_elem bits of any element
    let bitmask = (1 << bits_per_elem) - 1;

    // Write all the bits into the given bytestring
    let mut bit_idx = 0;
    while bit_idx < bits_per_elem * data.len() {
        let byte_idx = bit_idx / 8;
        let elem_idx = bit_idx / bits_per_elem;
        let bit_in_byte = bit_idx % 8;
        let bit_in_elem = bit_idx % bits_per_elem;

        // First clear the unused bits of the element we're going to write
        let elem_to_read = data[elem_idx] & bitmask;
        // Then  shift the element we're writing so the first unwritten bit is the lowest bit
        let elem_to_read = (elem_to_read >> bit_in_elem) as u8;

        // OR the bits into our byte, shifitng to align with the first unused bit in the byte
        out_buf[byte_idx] |= elem_to_read << bit_in_byte;

        // We just wrote either: 1) however many bits remained in the byte we were
        // reading, or 2) however many unused bits remained in the current element we were
        // writing to. Whichever is smaller.
        let just_wrote = core::cmp::min(8 - bit_in_byte, bits_per_elem - bit_in_elem);
        bit_idx += just_wrote
    }
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
        self.0.iter_mut().for_each(|w| *w &= MODULUS_Q - 1)
    }
}

impl<'a> Add for &'a RingElem {
    type Output = RingElem;

    fn add(self, other: &'a RingElem) -> Self::Output {
        poly_add(&self.0, &other.0)
    }
}

impl<'a> Sub for &'a RingElem {
    type Output = RingElem;

    fn sub(self, other: &'a RingElem) -> Self::Output {
        poly_sub(&self.0, &other.0)
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

            let prod_1 = &a * &b;
            let prod_2 = &b * &a;

            let sum_1 = &a + &b;
            let sum_2 = &b + &a;

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

            let kara_prod = karatsuba_mul_helper(&a.0, &b.0);
            let schoolbook_prod = &a * &b;

            assert_eq!(kara_prod, schoolbook_prod);
        }
    }

    #[test]
    fn from_bytes() {
        let mut rng = thread_rng();

        for _ in 0..1000 {
            // Check that from_bytes matches the reference impl from_bytes for N=2^13,2^10,2^1
            let bits_per_elem = 13;
            let mut bytes = vec![0u8; bits_per_elem * RING_DEG / 8];
            rng.fill_bytes(&mut bytes);
            assert_eq!(
                reference_impl_from_bytes_mod8192(&bytes),
                RingElem::from_bytes(&bytes, 13)
            );
            // Now check it matches the reference to_bytes impl
            let elem = RingElem::rand(&mut rng);
            let mut my_bytes = vec![0u8; bits_per_elem * RING_DEG / 8];
            let mut ref_bytes = vec![0u8; bits_per_elem * RING_DEG / 8];
            elem.to_bytes(&mut my_bytes, bits_per_elem);
            reference_impl_to_bytes_mod8192(&elem, &mut ref_bytes);
            assert_eq!(my_bytes, ref_bytes);

            let bits_per_elem = 10;
            let mut bytes = vec![0u8; bits_per_elem * RING_DEG / 8];
            rng.fill_bytes(&mut bytes);
            assert_eq!(
                reference_impl_from_bytes_mod1024(&bytes).0,
                RingElem::from_bytes(&bytes, 10).0,
            );

            let bits_per_elem = 1;
            let mut bytes = vec![0u8; bits_per_elem * RING_DEG / 8];
            rng.fill_bytes(&mut bytes);
            assert_eq!(
                reference_impl_from_bytes_mod2(&bytes).0,
                RingElem::from_bytes(&bytes, 1).0,
            );
            // Now check it matches the reference to_bytes impl
            let elem = RingElem::rand(&mut rng);
            let mut my_bytes = vec![0u8; bits_per_elem * RING_DEG / 8];
            let mut ref_bytes = vec![0u8; bits_per_elem * RING_DEG / 8];
            elem.to_bytes(&mut my_bytes, bits_per_elem);
            reference_impl_to_bytes_mod2(&elem, &mut ref_bytes);
            assert_eq!(my_bytes, ref_bytes);

            // Now check that to_bytes and from_bytes are inverses

            // Pick a random bits_per_elem
            let bits_per_elem = rng.gen_range(1..=13);
            let bitmask = (1 << bits_per_elem) - 1;

            // Generate a random element and make sure none of the values exceed 2^bits_per_elem
            let mut p = RingElem::rand(&mut rng);
            p.0.iter_mut().for_each(|e| *e &= bitmask);

            // Check that a round trip preserves the polynomial
            let mut p_bytes = vec![0u8; bits_per_elem * RING_DEG / 8];
            p.to_bytes(&mut p_bytes, bits_per_elem);
            assert_eq!(p, RingElem::from_bytes(&p_bytes, bits_per_elem));

            // Now other way around
            let mut p_bytes = vec![0u8; bits_per_elem * RING_DEG / 8];
            rng.fill_bytes(&mut p_bytes);
            let p = RingElem::from_bytes(&p_bytes, bits_per_elem);
            let mut new_p_bytes = vec![0u8; bits_per_elem * RING_DEG / 8];
            p.to_bytes(&mut new_p_bytes, bits_per_elem);
            assert_eq!(p_bytes, new_p_bytes);
        }
    }

    /// A nearly verbatim copy of the C reference impl of BS2POL_N where N = 2^13
    /// https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/pack_unpack.c#L101
    fn reference_impl_from_bytes_mod8192(b: &[u8]) -> RingElem {
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

        poly
    }

    /// A nearly verbatim copy of the C reference impl of BS2POL_N where N = 2^10
    /// https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/pack_unpack.c#L134
    fn reference_impl_from_bytes_mod1024(b: &[u8]) -> RingElem {
        let mut offset_byte;
        let mut offset_data;
        let mut poly = RingElem::default();
        let data = &mut poly.0;

        let bytes: Vec<u16> = b.iter().map(|&x| x as u16).collect();
        for j in 0..RING_DEG / 4 {
            offset_byte = 5 * j;
            offset_data = 4 * j;
            data[offset_data + 0] =
                (bytes[offset_byte + 0] & (0xff)) | ((bytes[offset_byte + 1] & 0x03) << 8);
            data[offset_data + 1] =
                ((bytes[offset_byte + 1] >> 2) & (0x3f)) | ((bytes[offset_byte + 2] & 0x0f) << 6);
            data[offset_data + 2] =
                ((bytes[offset_byte + 2] >> 4) & (0x0f)) | ((bytes[offset_byte + 3] & 0x3f) << 4);
            data[offset_data + 3] =
                ((bytes[offset_byte + 3] >> 6) & (0x03)) | ((bytes[offset_byte + 4] & 0xff) << 2);
        }

        poly
    }

    fn reference_impl_from_bytes_mod2(b: &[u8]) -> RingElem {
        let mut poly = RingElem::default();
        let data = &mut poly.0;

        let bytes: Vec<u16> = b.iter().map(|&x| x as u16).collect();
        for j in 0..32 {
            {
                for i in 0..8 {
                    data[j * 8 + i] = (bytes[j] >> i) & 0x01;
                }
            }
        }

        poly
    }

    fn reference_impl_to_bytes_mod8192(polyn: &RingElem, bytes: &mut [u8]) {
        let mut offset_byte: usize;
        let mut offset_data: usize;
        let data = polyn.0;

        for j in 0..RING_DEG / 8 {
            offset_byte = 13 * j;
            offset_data = 8 * j;
            bytes[offset_byte + 0] = (data[offset_data + 0] & (0xff)) as u8;
            bytes[offset_byte + 1] = ((data[offset_data + 0] >> 8) & 0x1f) as u8
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

    fn reference_impl_to_bytes_mod2(polyn: &RingElem, bytes: &mut [u8]) {
        let data = polyn.0;

        for j in 0..32 {
            for i in 0..8 {
                bytes[j] = bytes[j] | ((data[j * 8 + i] & 0x01) << i) as u8;
            }
        }
    }
}
