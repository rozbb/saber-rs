use crate::{
    consts::{MAX_MU, MODULUS_Q_BITS, RING_DEG},
    matrix_arith::Matrix,
    ring_arith::RingElem,
    util::deserialize,
};

use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128,
};

// Algorithm 16, GenSecret
/// Uses a random seed to generate an MLWR secret, i.e., an element in R^ℓ whose entries are
/// sampled according to a binomial distribution
pub(crate) fn gen_secret_from_seed<const L: usize, const MU: usize>(
    seed: &[u8; 32],
) -> Matrix<L, 1> {
    // Make a buffer of the correct length. We can't do const math there, so just make one of the
    // max length and then cut it down

    // Hash the seed and make an XOF
    let mut xof = {
        let mut h = Shake128::default();
        h.update(&seed[..]);
        h.finalize_xof()
    };

    // Output secret is L ring elements
    let mut ring_elems = [RingElem::default(); L];
    // Buffer to hold XOF bytes. Can't do const math here, so we make it the max size
    // and cut it down
    let mut backing_buf = [0u8; RING_DEG * MAX_MU / 8];
    let buf = &mut backing_buf[..RING_DEG * MU / 8];

    // Sample the secret
    for p in ring_elems.iter_mut() {
        // For each secret entry, read some XOF bytes into a buffer and parse it
        // as a ring elem whose coefficients are μ/2 bits
        xof.read(buf);
        // TODO optimization: deserialize is overkill here. we only need to count the
        // hamming weight of the deserialized words. This is done with the cheaper cbd()
        // function in the reference impl. It's pretty nasty-looking though
        let data: [u16; 2 * RING_DEG] = deserialize(buf, MU / 2);

        // Output element i is hamming(data[2*i]) - hamming(data[2*i+1])
        for (out_coeff, words) in p.0.iter_mut().zip(data.chunks_exact(2)) {
            let hamming1 = words[0].count_ones() as u16;
            let hamming2 = words[1].count_ones() as u16;
            *out_coeff = hamming1.wrapping_sub(hamming2);
        }
    }

    Matrix([ring_elems]).transpose()
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
