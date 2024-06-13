use crate::{
    matrix_arith::Matrix,
    ring_arith::{deserialize, RingElem, MODULUS_Q_BITS, RING_DEG},
};

use rand_core::CryptoRngCore;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128,
};

/// The maximum possible μ value is μ=10, set by Lightsaber
const MAX_MU: usize = 10;

/// The maximum possible l value is l=4, set by FireSaber
const MAX_L: usize = 4;

// Algorithm 16, GenSecret
/// Uses a random seed to generate an MLWR secret, i.e., an element in R^ℓ whose entries are
/// sampled according to a binomial distribution
pub(crate) fn gen_secret_from_seed<const L: usize, const MU: usize>(
    seed: &[u8; 32],
) -> Matrix<L, 1> {
    // Make a buffer of the correct length. We can't do const math there, so just make one of the
    // max length and then cut it down
    let mut backing_buf = [0u8; MAX_L * RING_DEG * MAX_MU / 8];
    let buf = &mut backing_buf[..L * RING_DEG * MU / 8];

    // Hash seed to fill up buf with randomness
    Shake128::digest_xof(&seed, buf);

    // Create the L secret polynomials
    let mut polyns = [RingElem::default(); L];
    for (polyn_idx, chunk) in buf.chunks(RING_DEG * MU / 8).enumerate() {
        let data: [u16; 2 * RING_DEG] = deserialize(chunk, MU / 2);
        // Our output element is hamming(data[i]) - hamming(data[i+1]) for every i, skipping by 2
        let polyn = &mut polyns[polyn_idx];
        for i in 0..RING_DEG {
            let hamming1 = data[2 * i].count_ones() as u16;
            let hamming2 = data[2 * i + 1].count_ones() as u16;
            polyn.0[i] = hamming1.wrapping_sub(hamming2);
        }
        // Mask off the high order bits
        polyn.reduce_mod_2pow(MODULUS_Q_BITS);
    }

    Matrix([polyns])
}

// Algorithm 15, GenMatrix
/// Uses a random seed to generate a uniform matrix in R^{ℓ×ℓ}
pub(crate) fn gen_matrix_from_seed<const L: usize>(seed: &[u8; 32]) -> Matrix<L, L> {
    // Make a buffer of the correct length. We can't do const math there, so just make one of the
    // max length and then cut it down
    let mut backing_buf = [0u8; MAX_L * MAX_L * RING_DEG * MODULUS_Q_BITS / 8];
    let buf = &mut backing_buf[..L * L * RING_DEG * MODULUS_Q_BITS / 8];

    // Hash seed to fill up buf with randomness
    Shake128::digest_xof(&seed, buf);

    let mut mat = Matrix::default();
    for (idx, chunk) in buf.chunks(RING_DEG * MODULUS_Q_BITS / 8).enumerate() {
        let polyn = RingElem::from_bytes(chunk, MODULUS_Q_BITS);

        let i = idx / L;
        let j = idx % L;
        mat.0[i][j] = polyn;
    }

    mat
}

#[test]
fn test_gen_secret() {
    use rand::RngCore;
    let mut rng = rand::thread_rng();
    const L: usize = 4;
    const MU: usize = 10;
    let mut seed = [0u8; 32];
    rng.fill_bytes(&mut seed);
    let out = gen_secret_from_seed::<L, MU>(&seed);
}
