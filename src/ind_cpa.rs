use crate::{
    gen::{gen_matrix_from_seed, gen_secret_from_seed},
    matrix_arith::Matrix,
    ring_arith::{deserialize, serialize, RingElem, MODULUS_P_BITS, MODULUS_Q_BITS, RING_DEG},
};

use rand_core::CryptoRngCore;

const H1_VAL: u16 = 1 << (MODULUS_Q_BITS - MODULUS_P_BITS - 1);

struct IndCpaSecretKey<const L: usize>(Matrix<L, 1>);
struct IndCpaPublicKey<const L: usize> {
    seed: [u8; 32],
    vec: Matrix<L, 1>,
}

impl<const L: usize> IndCpaSecretKey<L> {
    fn serialize(&self, out_buf: &mut [u8]) {
        assert_eq!(out_buf.len(), L * MODULUS_Q_BITS * RING_DEG / 8);
        self.0.to_bytes(out_buf, MODULUS_Q_BITS);
    }
}

impl<const L: usize> IndCpaPublicKey<L> {
    fn serialize(&self, out_buf: &mut [u8]) {
        assert_eq!(out_buf.len(), 32 + L * MODULUS_P_BITS * RING_DEG / 8);

        // Write out the pubkey seed
        out_buf[..32].copy_from_slice(&self.seed);
        // Write out the rest of the pubkey
        self.vec.to_bytes(&mut out_buf[32..], MODULUS_P_BITS);
    }
}

// Algorithm 17, Saber.PKE.KeyGen
/// Generates a keypair with a secret from R^ℓ with bionimal parameter μ
fn gen_keypair<const L: usize, const MU: usize>(
    rng: &mut impl CryptoRngCore,
) -> (IndCpaSecretKey<L>, IndCpaPublicKey<L>) {
    let mut matrix_seed = [0u8; 32];
    let mut secret_seed = [0u8; 32];
    rng.fill_bytes(&mut matrix_seed);
    rng.fill_bytes(&mut secret_seed);

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

fn dec<const L: usize, const MODULUS_T_BITS: usize>(
    sk: &IndCpaSecretKey<L>,
    ciphertext: &[u8],
) -> [u8; 32] {
    let (c_bytes, bprime_bytes) = ciphertext.split_at(MODULUS_T_BITS * RING_DEG / 8);

    let bprime: Matrix<L, 1> = Matrix::from_bytes(bprime_bytes, MODULUS_P_BITS);

    let mut c = RingElem::from_bytes(c_bytes, MODULUS_T_BITS);
    c.shift_left(MODULUS_P_BITS - MODULUS_T_BITS);

    let v = bprime.mul_transpose(&sk.0);
    let v = v.0[0][0];

    // Compute v - c + h₂
    let mut c = &v - &c;
    let h2_val = (1 << (MODULUS_P_BITS - 2)) - (1 << (MODULUS_P_BITS - MODULUS_T_BITS - 1))
        + (1 << (MODULUS_Q_BITS - MODULUS_P_BITS - 1));
    c.wrapping_add_to_all(h2_val);
    c.shift_right(MODULUS_P_BITS - 1);

    let mut m = [0u8; 32];
    c.to_bytes(&mut m, 1);
    m
}

fn enc_deterministic<const L: usize, const MU: usize, const MODULUS_T_BITS: usize>(
    pk: &IndCpaPublicKey<L>,
    sk: &IndCpaSecretKey<L>,
    msg: &[u8; 32],
    seed: &[u8; 32],
    out_buf: &mut [u8],
) {
    assert_eq!(
        out_buf.len(),
        MODULUS_T_BITS * RING_DEG / 8 + L * MODULUS_P_BITS * RING_DEG / 8
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

    c.to_bytes(
        &mut out_buf[..MODULUS_T_BITS * RING_DEG / 8],
        MODULUS_T_BITS,
    );
    bprime.to_bytes(
        &mut out_buf[MODULUS_T_BITS * RING_DEG / 8..],
        MODULUS_P_BITS,
    );
}

#[cfg(test)]
mod test {
    use super::*;

    use rand::RngCore;

    // Tests that Dec(Enc(m)) == m
    #[test]
    fn enc_dec() {
        const L: usize = 2;
        const MODULUS_T_BITS: usize = 3;
        const MU: usize = 10;

        let mut rng = rand::thread_rng();

        let (sk, pk) = gen_keypair::<L, MU>(&mut rng);

        let mut enc_seed = [0u8; 32];
        let mut msg = [0u8; 32];
        rng.fill_bytes(&mut enc_seed);
        rng.fill_bytes(&mut msg);
        let mut ct_buf = [0u8; MODULUS_T_BITS * RING_DEG / 8 + L * MODULUS_P_BITS * RING_DEG / 8];

        enc_deterministic::<L, MU, MODULUS_T_BITS>(&pk, &sk, &msg, &enc_seed, &mut ct_buf);
        let recovered_msg = dec::<L, MODULUS_T_BITS>(&sk, &ct_buf);
        assert_eq!(msg, recovered_msg);
    }
}
