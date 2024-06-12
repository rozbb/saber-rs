use crate::{
    gen::{gen_matrix_from_seed, gen_secret_from_seed},
    matrix_arith::Matrix,
    ring_arith::{MODULUS_P_BITS, MODULUS_Q_BITS, RING_DEG},
};

use rand_core::CryptoRngCore;

// Algorithm 17, Saber.PKE.KeyGen
/// Generates a keypair with a secret from R^ℓ with bionimal parameter μ
fn gen_keypair<const L: usize, const MU: usize>(
    rng: &mut impl CryptoRngCore,
    sk_out: &mut [u8],
    pk_out: &mut [u8],
) {
    assert_eq!(sk_out.len(), L * MODULUS_Q_BITS * RING_DEG / 8);
    assert_eq!(pk_out.len(), 32 + L * MODULUS_P_BITS * RING_DEG / 8);

    let mut matrix_seed = [0u8; 32];
    let mut secret_seed = [0u8; 32];
    rng.fill_bytes(&mut matrix_seed);
    rng.fill_bytes(&mut secret_seed);

    let mat_a = gen_matrix_from_seed::<L>(matrix_seed);
    let vec_s = gen_secret_from_seed::<L, MU>(secret_seed);
    let b = {
        let mut prod = mat_a.mul_transpose(&vec_s);
        // Add the h vector, which consists enitrely of 4's
        prod.wrapping_add_to_all(4);
        // Now shift all the coefficients by EQ - EP = 13 - 10 = 3
        prod.shift_right(3);
        prod
    };

    // Write out the secret
    vec_s.to_bytes(sk_out, MODULUS_Q_BITS);

    // Write out the pubkey seed
    pk_out[..32].copy_from_slice(&matrix_seed);
    // Write out the rest of the pubkey
    b.to_bytes(&mut pk_out[32..], MODULUS_P_BITS);
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn keygen() {
        const L: usize = 4;
        const MU: usize = 10;

        let mut rng = rand::thread_rng();

        let mut sk_buf = [0u8; L * MODULUS_Q_BITS * RING_DEG / 8];
        let mut pk_buf = [0u8; 32 + L * MODULUS_P_BITS * RING_DEG / 8];
        gen_keypair::<L, MU>(&mut rng, &mut sk_buf, &mut pk_buf);
    }
}
