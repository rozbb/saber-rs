use crate::{
    gen::{gen_matrix_from_seed, gen_secret_from_seed},
    matrix_arith::Matrix,
    ring_arith::{deserialize, serialize, RingElem, MODULUS_P_BITS, MODULUS_Q_BITS, RING_DEG},
};

use rand_core::CryptoRngCore;
use sha3::{digest::ExtendableOutput, Shake128};

const H1_VAL: u16 = 1 << (MODULUS_Q_BITS - MODULUS_P_BITS - 1);

pub(crate) struct IndCpaSecretKey<const L: usize>(Matrix<L, 1>);

#[derive(Clone)]
pub struct IndCpaPublicKey<const L: usize> {
    seed: [u8; 32],
    vec: Matrix<L, 1>,
}

impl<const L: usize> IndCpaSecretKey<L> {
    pub(crate) fn to_bytes(&self, out_buf: &mut [u8]) {
        assert_eq!(out_buf.len(), Self::serialized_len());

        self.0.to_bytes(out_buf, MODULUS_Q_BITS);
    }

    pub(crate) fn serialized_len() -> usize {
        L * MODULUS_Q_BITS * RING_DEG / 8
    }
}

impl<const L: usize> IndCpaPublicKey<L> {
    pub(crate) fn to_bytes(&self, out_buf: &mut [u8]) {
        let out_size = Self::serialized_len();
        assert_eq!(out_buf.len(), out_size);

        // We write out the LWR sample and then the seed. The spec actually said to do the opposite, but the reference impl does it this way
        // https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/SABER_indcpa.c#L42

        // Write out the LWR sample
        self.vec
            .to_bytes(&mut out_buf[..out_size - 32], MODULUS_P_BITS);
        // Write out the pubkey seed
        out_buf[out_size - 32..].copy_from_slice(&self.seed);
    }

    pub(crate) fn serialized_len() -> usize {
        32 + L * MODULUS_P_BITS * RING_DEG / 8
    }
}

// Algorithm 17, Saber.PKE.KeyGen
/// Generates a keypair with a secret from R^ℓ with bionimal parameter μ
pub(crate) fn gen_keypair<const L: usize, const MU: usize>(
    rng: &mut impl CryptoRngCore,
) -> (IndCpaSecretKey<L>, IndCpaPublicKey<L>) {
    let mut matrix_seed = [0u8; 32];
    let mut matrix_seed_unhashed = [0u8; 32];
    let mut secret_seed = [0u8; 32];
    rng.fill_bytes(&mut matrix_seed_unhashed);
    rng.fill_bytes(&mut secret_seed);

    // Before using the matrix seed, we need to hash it
    Shake128::digest_xof(&matrix_seed_unhashed, &mut matrix_seed);

    let mat_a = gen_matrix_from_seed::<L>(&matrix_seed);
    let vec_s = gen_secret_from_seed::<L, MU>(&secret_seed);
    let b = {
        let mut prod = mat_a.mul_transpose(&vec_s);
        prod.wrapping_add_to_all(H1_VAL);
        prod.shift_right(MODULUS_Q_BITS - MODULUS_P_BITS);
        prod
    };

    (
        IndCpaSecretKey(vec_s),
        IndCpaPublicKey {
            seed: matrix_seed,
            vec: b,
        },
    )
}

pub(crate) fn dec<const L: usize, const MODULUS_T_BITS: usize>(
    sk: &IndCpaSecretKey<L>,
    ciphertext: &[u8],
) -> [u8; 32] {
    assert_eq!(
        ciphertext.len(),
        L * MODULUS_P_BITS * RING_DEG / 8 + MODULUS_T_BITS * RING_DEG / 8
    );
    let (bprime_bytes, c_bytes) = ciphertext.split_at(L * MODULUS_P_BITS * RING_DEG / 8);

    let bprime: Matrix<L, 1> = Matrix::from_bytes(bprime_bytes, MODULUS_P_BITS);

    let mut c = RingElem::from_bytes(c_bytes, MODULUS_T_BITS);
    c.shift_left(MODULUS_P_BITS - MODULUS_T_BITS);

    let v = bprime.mul_transpose(&sk.0);
    let v = v.0[0][0];

    // Compute v - c + h₂
    let mut mprime = &v - &c;
    let h2_val = (1 << (MODULUS_P_BITS - 2)) - (1 << (MODULUS_P_BITS - MODULUS_T_BITS - 1))
        + (1 << (MODULUS_Q_BITS - MODULUS_P_BITS - 1));
    mprime.wrapping_add_to_all(h2_val);
    mprime.shift_right(MODULUS_P_BITS - 1);

    let mut m = [0u8; 32];
    mprime.to_bytes(&mut m, 1);
    m
}

pub(crate) const fn ciphertext_size<
    const L: usize,
    const MU: usize,
    const MODULUS_T_BITS: usize,
>() -> usize {
    MODULUS_T_BITS * RING_DEG / 8 + L * MODULUS_P_BITS * RING_DEG / 8
}

pub(crate) fn enc_deterministic<const L: usize, const MU: usize, const MODULUS_T_BITS: usize>(
    pk: &IndCpaPublicKey<L>,
    msg: &[u8; 32],
    seed: &[u8; 32],
    out_buf: &mut [u8],
) {
    assert_eq!(
        out_buf.len(),
        L * MODULUS_P_BITS * RING_DEG / 8 + MODULUS_T_BITS * RING_DEG / 8
    );

    let mat_a = gen_matrix_from_seed::<L>(&pk.seed);
    let vec_sprime = gen_secret_from_seed::<L, MU>(seed);

    let bprime = {
        let mut prod = mat_a.mul(&vec_sprime);
        prod.wrapping_add_to_all(H1_VAL);
        prod.shift_right(MODULUS_Q_BITS - MODULUS_P_BITS);
        prod
    };

    let vprime: Matrix<1, 1> = pk.vec.mul_transpose(&vec_sprime);
    let vprime = vprime.0[0][0];

    let mut msg_polyn = RingElem(deserialize(msg, 1));
    msg_polyn.shift_left(MODULUS_P_BITS - 1);

    // Compute v' - mp + h₁
    let mut c = &vprime - &msg_polyn;
    c.wrapping_add_to_all(H1_VAL);
    c.shift_right(MODULUS_P_BITS - MODULUS_T_BITS);

    // Serialization order accoring to the reference implementation is b' || c
    // https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/SABER_indcpa.c#L68-L79
    let (bprime_buf, c_buf) = out_buf.split_at_mut(L * MODULUS_P_BITS * RING_DEG / 8);
    bprime.to_bytes(bprime_buf, MODULUS_P_BITS);
    c.to_bytes(c_buf, MODULUS_T_BITS);
}

#[cfg(test)]
mod test {
    use super::*;

    use rand::RngCore;

    // Tests that Dec(Enc(m)) == m
    #[test]
    fn enc_dec() {
        test_enc_dec::<2, 3, 10>();
        test_enc_dec::<3, 4, 8>();
        test_enc_dec::<4, 6, 6>();
    }

    fn test_enc_dec<const L: usize, const MODULUS_T_BITS: usize, const MU: usize>() {
        let mut rng = rand::thread_rng();

        for _ in 0..100 {
            let (sk, pk) = gen_keypair::<L, MU>(&mut rng);

            let mut enc_seed = [0u8; 32];
            let mut msg = [0u8; 32];
            rng.fill_bytes(&mut enc_seed);
            rng.fill_bytes(&mut msg);
            let mut ct_buf =
                vec![0u8; MODULUS_T_BITS * RING_DEG / 8 + L * MODULUS_P_BITS * RING_DEG / 8];

            enc_deterministic::<L, MU, MODULUS_T_BITS>(&pk, &msg, &enc_seed, &mut ct_buf);
            let recovered_msg = dec::<L, MODULUS_T_BITS>(&sk, &ct_buf);
            assert_eq!(msg, recovered_msg);
        }
    }
}
