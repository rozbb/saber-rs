//! This file implements methods for generating uniform matrices and binomially distributed vectors

use crate::{
    arithmetic::{Matrix, RingElem},
    consts::{DOMSEP_GENMAT, DOMSEP_GENSEC, MAX_MU, MODULUS_Q_BITS, RING_DEG},
};

use turboshake::digest::{ExtendableOutput, Update, XofReader};
use turboshake::{CTurboShake128, CTurboShake256};

/// Computes the Centered Binomial Distribution using the given bytes as randomness
fn cbd<const MU: usize>(buf: &[u8], out: &mut RingElem) {
    let half = MU / 2;
    let mask: u32 = (1 << half) - 1;

    // Special case for Kopis-768: Each coefficient uses exactly 1 byte, with the low nibble
    // as the positive half and the high nibble as the negative half. So we don't need a
    // buffer to read bits
    if MU == 8 {
        // Specialized fast path for Kopis-768 (MU=8): each coefficient = one byte, no bit shifting
        for (coeff, &byte) in out.0.iter_mut().zip(buf.iter()) {
            let a = (byte & 0x0F).count_ones() as u16;
            let b = (byte >> 4).count_ones() as u16;
            *coeff = a.wrapping_sub(b);
        }
    } else {
        // General path for MU=6 (Kopis-1024) and MU=10 (Kopis-512).
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

/// Uses a random seed to generate an MLWR secret, i.e., an element in R^ℓ whose entries are
/// sampled according to a binomial distribution.
pub(crate) fn gen_secret_from_seed<const L: usize, const MU: usize>(
    seed: &[u8; 32],
) -> Matrix<L, 1> {
    let mut secret = Matrix::default();
    // Buffer to hold XOF bytes. Can't do const math here, so we make it the max size
    // and cut it down
    let mut backing_buf = [0u8; RING_DEG * MAX_MU / 8];
    let buf = &mut backing_buf[..RING_DEG * MU / 8];

    // Sample the secret using the Centered Binomial Distribution
    for (i, row) in secret.0.iter_mut().enumerate() {
        let mut hasher = CTurboShake256::<DOMSEP_GENSEC>::default();
        hasher.update(seed);
        hasher.update(&[i as u8]);
        let mut reader = hasher.finalize_xof();
        reader.read(buf);
        cbd::<MU>(buf, &mut row[0]);
    }

    secret
}

/// Uses a random seed to generate a uniform matrix in R^{ℓ×ℓ}.
///
/// For each element (i,j), we compute TurboSHAKE128(seed || i || j, 256*13/8, DOMSEP_GENMAT).
pub(crate) fn gen_matrix_from_seed<const L: usize>(seed: &[u8; 32]) -> Matrix<L, L> {
    // Our output is a matrix of ring elements
    let mut mat = Matrix::default();
    // For each ring element we need to sample the same number of bytes
    let mut buf = [0u8; RING_DEG * MODULUS_Q_BITS / 8];

    // Construct the matrix entries
    for (i, row) in mat.0.iter_mut().enumerate() {
        for (j, p) in row.iter_mut().enumerate() {
            let mut hasher = CTurboShake128::<DOMSEP_GENMAT>::default();
            hasher.update(seed);
            hasher.update(&[i as u8]);
            hasher.update(&[j as u8]);
            let mut reader = hasher.finalize_xof();
            reader.read(&mut buf);
            *p = RingElem::from_bytes(&buf, MODULUS_Q_BITS);
        }
    }

    mat
}
