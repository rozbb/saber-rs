//! This file implements methods for generating uniform matrices and binomially distributed vectors

use crate::{
    arithmetic::{Matrix, RingElem},
    consts::{MAX_MU, MODULUS_Q_BITS, RING_DEG},
};

use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128,
};

/// Computes the Centered Binomial Distribution directly from raw bytes.
///
/// Each ring coefficient is sampled from CBD(μ/2): take μ random bits, split into two halves
/// of μ/2 bits each, and output popcount(first_half) - popcount(second_half).
///
/// For MU=8 (Saber), this is byte-aligned: each coefficient uses exactly 1 byte, with the
/// low nibble as the positive half and the high nibble as the negative half.
///
/// For MU=6 (FireSaber) and MU=10 (LightSaber), we read MU bits at a time from the byte
/// stream using bitwise extraction.
fn cbd<const MU: usize>(buf: &[u8], out: &mut RingElem) {
    let half = MU / 2;
    let mask: u32 = (1 << half) - 1;

    if MU == 8 {
        // Specialized fast path for Saber (MU=8): each coefficient = one byte, no bit shifting
        for (coeff, &byte) in out.0.iter_mut().zip(buf.iter()) {
            let a = (byte & 0x0F).count_ones() as u16;
            let b = (byte >> 4).count_ones() as u16;
            *coeff = a.wrapping_sub(b);
        }
    } else {
        // General path for MU=6 (FireSaber) and MU=10 (LightSaber).
        // Read MU bits at a time, spanning up to 3 bytes when not byte-aligned.
        let mut bit_pos = 0;
        for coeff in out.0.iter_mut() {
            let byte_idx = bit_pos / 8;
            let bit_in_byte = bit_pos % 8;

            // Read up to 3 bytes to cover MU bits starting at bit_in_byte.
            // Worst case: MU=10 starting at bit 7 needs bits 7..16, spanning 3 bytes.
            let mut raw: u32 = buf[byte_idx] as u32;
            if byte_idx + 1 < buf.len() {
                raw |= (buf[byte_idx + 1] as u32) << 8;
            }
            if byte_idx + 2 < buf.len() {
                raw |= (buf[byte_idx + 2] as u32) << 16;
            }
            raw >>= bit_in_byte;

            let a = (raw & mask).count_ones() as u16;
            let b = ((raw >> half) & mask).count_ones() as u16;
            *coeff = a.wrapping_sub(b);

            bit_pos += MU;
        }
    }
}

// Algorithm 16, GenSecret
/// Uses a random seed to generate an MLWR secret, i.e., an element in R^ℓ whose entries are
/// sampled according to a binomial distribution
pub(crate) fn gen_secret_from_seed<const L: usize, const MU: usize>(
    seed: &[u8; 32],
) -> Matrix<L, 1> {
    // Hash the seed and make an XOF
    let mut xof = {
        let mut h = Shake128::default();
        h.update(&seed[..]);
        h.finalize_xof()
    };

    // Output secret is a column vector (Matrix<L, 1>). We construct it directly
    // rather than building a row vector and transposing, to avoid copying L RingElems.
    let mut secret = Matrix::default();
    // Buffer to hold XOF bytes. Can't do const math here, so we make it the max size
    // and cut it down
    let mut backing_buf = [0u8; RING_DEG * MAX_MU / 8];
    let buf = &mut backing_buf[..RING_DEG * MU / 8];

    // Sample the secret using the Centered Binomial Distribution
    for row in secret.0.iter_mut() {
        xof.read(buf);
        cbd::<MU>(buf, &mut row[0]);
    }

    secret
}

// Algorithm 15, GenMatrix
/// Uses a random seed to generate a uniform matrix in R^{ℓ×ℓ}
pub(crate) fn gen_matrix_from_seed<const L: usize>(seed: &[u8; 32]) -> Matrix<L, L> {
    // Hash the seed and make an XOF
    let mut xof = {
        let mut h = Shake128::default();
        h.update(&seed[..]);
        h.finalize_xof()
    };

    // Our output is a matrix of ring elements
    let mut mat = Matrix::default();
    // For each ring element we need to sample the same number of bytes
    let mut buf = [0u8; RING_DEG * MODULUS_Q_BITS / 8];

    // Construct the matrix entries
    for row in mat.0.iter_mut() {
        for p in row.iter_mut() {
            // For each matrix entry, read some XOF bytes into a buffer and parse it
            // as a ring elem whose coefficients are all q bits
            xof.read(&mut buf);
            *p = RingElem::from_bytes(&buf, MODULUS_Q_BITS);
        }
    }

    mat
}
