use core::ops::Mul;

use rand_core::CryptoRngCore;

/// The modulus of our base ring Z/2^13 Z
const MODULUS: u16 = 1 << 13;

/// The value of -1 in Z/2^13 Z
const NEG_ONE: u16 = MODULUS - 1;

/// The degree of the polynomial ring over Z/2^13 Z
const RING_DEG: usize = 256;

/// An element of the ring (Z/2^13 Z)[X] / (X^256 + 1)
// INVARIANT: The top 3 bits of each coefficient is always 0
#[derive(Eq, PartialEq, Debug, Clone)]
struct RingElem([u16; RING_DEG]);

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
                let (idx, idx_is_past_256) = (i as u8).overflowing_add(j as u8);
                // The coeff a·X^{256 + i} equals -a·X^i since X^256 + 1 == 0 in our polyn ring
                if idx_is_past_256 {
                    prod = prod.wrapping_mul(NEG_ONE);
                }
                result[idx as usize] = result[idx as usize].wrapping_add(prod);
            }
        }

        // Now reduce everything
        for i in 0..RING_DEG {
            result[i] %= MODULUS;
        }

        RingElem(result)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use rand::thread_rng;

    // Checks that a * b == b * a for ring elements a, b
    #[test]
    fn mult_commutativity() {
        let mut rng = thread_rng();

        for _ in 0..100 {
            let a = RingElem::rand(&mut rng);
            let b = RingElem::rand(&mut rng);

            assert_eq!(&a * &b, &b * &a);
        }
    }
}
