//! This file implements the IND-CPA-secure Kopis PKE scheme

use crate::{
    arithmetic::{Matrix, RingElem},
    consts::{
        DOMSEP_KGEXPAND, DOMSEP_PKHASH, MAX_L, MAX_T, MODULUS_P_BITS, MODULUS_Q_BITS, RING_DEG,
    },
    gen::{gen_matrix_from_seed, gen_secret_from_seed},
    ser::deserialize_generic,
    turboshake256_hash,
};

use turboshake::digest::{ExtendableOutput, Update, XofReader};
use turboshake::CTurboShake256;

const H1_VAL: u16 = 1 << (MODULUS_Q_BITS - MODULUS_P_BITS - 1);

/// A secret key for the IND-CPA-secure Kopis PKE scheme (expanded form)
pub(crate) struct PkeSecretKey<const L: usize>(Matrix<L, 1>);

/// A public key for the IND-CPA-secure Kopis PKE scheme
#[derive(Clone)]
pub struct PkePublicKey<const L: usize> {
    matrix_seed: [u8; 32],
    vec: Matrix<L, 1>,
    /// The public matrix A, expanded from `matrix_seed`. Cached here so that repeated
    /// encryptions (e.g. every encapsulation and every FO re-encryption during decapsulation)
    /// don't have to re-run the XOF that derives A from the seed. This mirrors the "unpacked"
    /// public-key form used by other KEM implementations. It is never serialized: `serialize`
    /// still writes only `vec || matrix_seed`, and `from_bytes` re-derives it.
    mat_a: Matrix<L, L>,
}

impl<const L: usize> PkePublicKey<L> {
    pub const SERIALIZED_LEN: usize = 32 + L * MODULUS_P_BITS * RING_DEG / 8;

    /// Serializes this public key to a byte string. `out_buf` MUST have length SERIALIZED_LEN
    pub(crate) fn serialize(&self, out_buf: &mut [u8]) {
        let out_size = Self::SERIALIZED_LEN;
        assert_eq!(out_buf.len(), out_size);

        // Write out the LWR sample, then the seed
        self.vec
            .serialize(&mut out_buf[..out_size - 32], MODULUS_P_BITS);
        // Write out the pubkey seed
        out_buf[out_size - 32..].copy_from_slice(&self.matrix_seed);
    }

    pub(crate) fn from_bytes(bytes: &[u8]) -> Self {
        assert_eq!(bytes.len(), Self::SERIALIZED_LEN);

        let (vec_bytes, seed) = bytes.split_at(Self::SERIALIZED_LEN - 32);
        let vec = Matrix::deserialize_10(vec_bytes);
        let matrix_seed: [u8; 32] = seed.try_into().unwrap(); // checked above
        let mat_a = gen_matrix_from_seed::<L>(&matrix_seed);
        Self {
            matrix_seed,
            vec,
            mat_a,
        }
    }

    /// Returns the public key hash
    pub(crate) fn hash(&self) -> [u8; 32] {
        // pkh = TurboSHAKE256(pk, 32, DOMSEP_PKHASH)
        let mut buf = [0u8; max_pke_pubkey_serialized_len()];
        let pk_slice = &mut buf[..PkePublicKey::<L>::SERIALIZED_LEN];
        self.serialize(pk_slice);
        turboshake256_hash::<DOMSEP_PKHASH>(pk_slice, &[])
    }
}

/// The maximum length of a serialized public key, for all parameter choices
pub(crate) const fn max_pke_pubkey_serialized_len() -> usize {
    32 + MAX_L * MODULUS_P_BITS * RING_DEG / 8
}

/// The maximum length of a ciphertext (PKE or KEM, since they're the same), for all parameter
/// choices, for a message that is 32-bytes.
pub const fn max_ciphertext_len() -> usize {
    // b' is in R^l_P and c is in R_T
    MAX_T * RING_DEG / 8 + MAX_L * MODULUS_P_BITS * RING_DEG / 8
}

/// The length of a ciphertext (PKE or KEM, since they're the same) for a given parameter choice,
/// for a message that is 32-bytes.
pub const fn ciphertext_len<const L: usize, const T: usize>() -> usize {
    // b' is in R^l_P and c is in R_T
    L * MODULUS_P_BITS * RING_DEG / 8 + T * RING_DEG / 8
}

/// Expands a 32-byte secret key into the full decapsulation key components.
///
/// Returns (vec_s, z, pk, pkh) where:
/// - vec_s is the secret vector (as PkeSecretKey)
/// - z is 32 bytes used for rejection in decapsulation
/// - pk is the public key
/// - pkh is the hash of the public key
pub(crate) fn expand_decap_key<const L: usize, const MU: usize>(
    sk: &[u8; 32],
) -> (PkeSecretKey<L>, [u8; 32], PkePublicKey<L>, [u8; 32]) {
    // mat_seed || secret_seed || r = TurboSHAKE256(sk || L, 96, DOMSEP_KGEXPAND)
    let mut mat_seed = [0u8; 32];
    let mut secret_seed = [0u8; 32];
    let mut z = [0u8; 32];

    let mut xof = {
        let mut hasher = CTurboShake256::<DOMSEP_KGEXPAND>::default();
        hasher.update(sk);
        hasher.update(&[L as u8]);
        hasher.finalize_xof()
    };
    xof.read(&mut mat_seed);
    xof.read(&mut secret_seed);
    xof.read(&mut z);

    let mat_a = gen_matrix_from_seed::<L>(&mat_seed);
    let vec_s = gen_secret_from_seed::<L, MU>(&secret_seed);

    // vec_b = RoundToR10(transpose(mat_A) * vec_s)
    let b = {
        let mut prod = mat_a.mul_transpose(&vec_s);
        prod.wrapping_add_to_all(H1_VAL);
        prod.shift_right(MODULUS_Q_BITS - MODULUS_P_BITS);
        prod
    };

    let pk = PkePublicKey {
        matrix_seed: mat_seed,
        vec: b,
        mat_a,
    };
    let pkh = pk.hash();

    (PkeSecretKey(vec_s), z, pk, pkh)
}

/// Decrypts a ciphertext using the given secret key. `ciphertext` MUST have length
/// `ciphertext_len::<L, T>()`.
pub(crate) fn decrypt<const L: usize, const T: usize>(
    sk: &PkeSecretKey<L>,
    ciphertext: &[u8],
) -> [u8; 32] {
    assert_eq!(ciphertext.len(), ciphertext_len::<L, T>());
    // b' is in R^l_P and c is in R_T
    let (bprime_bytes, c_bytes) = ciphertext.split_at(L * MODULUS_P_BITS * RING_DEG / 8);

    let bprime: Matrix<L, 1> = Matrix::deserialize_10(bprime_bytes);

    let mut c = RingElem::deserialize(c_bytes, T);
    c.shift_left(MODULUS_P_BITS - T);

    let v = bprime.mul_transpose(&sk.0);
    let v = v.0[0][0];

    // Compute v - c + h₂
    let mut mprime = &v - &c;
    let h2_val = (1 << (MODULUS_P_BITS - 2)) - (1 << (MODULUS_P_BITS - T - 1))
        + (1 << (MODULUS_Q_BITS - MODULUS_P_BITS - 1));
    mprime.wrapping_add_to_all(h2_val);
    mprime.shift_right(MODULUS_P_BITS - 1);

    let mut m = [0u8; 32];
    mprime.serialize(&mut m, 1);
    m
}

/// Encrypts a message with a given public key and randomness (`coins`).
/// `out_buf` MUST have length `ciphertext_len::<L, T>()`.
pub(crate) fn encrypt_deterministic<const L: usize, const MU: usize, const T: usize>(
    pk: &PkePublicKey<L>,
    msg: &[u8; 32],
    coins: &[u8; 32],
    out_buf: &mut [u8],
) {
    assert_eq!(out_buf.len(), ciphertext_len::<L, T>());

    let vec_sprime = gen_secret_from_seed::<L, MU>(coins);

    let bprime = {
        let mut prod = pk.mat_a.mul(&vec_sprime);
        prod.wrapping_add_to_all(H1_VAL);
        prod.shift_right(MODULUS_Q_BITS - MODULUS_P_BITS);
        prod
    };

    let vprime: Matrix<1, 1> = pk.vec.mul_transpose(&vec_sprime);
    let vprime = vprime.0[0][0];

    let mut msg_polyn = RingElem(deserialize_generic(msg, 1));
    msg_polyn.shift_left(MODULUS_P_BITS - 1);

    // Compute v' - mp + h₁
    let mut c = &vprime - &msg_polyn;
    c.wrapping_add_to_all(H1_VAL);
    c.shift_right(MODULUS_P_BITS - T);

    // b' is in R^l_P and c is in R_T
    let (bprime_buf, c_buf) = out_buf.split_at_mut(L * MODULUS_P_BITS * RING_DEG / 8);
    bprime.serialize(bprime_buf, MODULUS_P_BITS);
    c.serialize(c_buf, T);
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::consts::*;

    use rand::RngCore;

    // Tests that Decrypt(Encrypt(m)) == m
    #[test]
    fn encryption_correctness() {
        test_enc_dec::<KOPIS512_L, KOPIS512_T, KOPIS512_MU>();
        test_enc_dec::<KOPIS768_L, KOPIS768_T, KOPIS768_MU>();
        test_enc_dec::<KOPIS1024_L, KOPIS1024_T, KOPIS1024_MU>();
    }

    // Helper function that encrypts and decrypts a random 32-byte message
    fn test_enc_dec<const L: usize, const T: usize, const MU: usize>() {
        let mut rng = rand::rng();
        let mut backing_buf = [0u8; MAX_T * RING_DEG / 8 + MAX_L * MODULUS_P_BITS * RING_DEG / 8];

        for _ in 0..100 {
            // Generate a random secret key seed and expand it
            let mut sk_seed = [0u8; 32];
            rng.fill_bytes(&mut sk_seed);
            let (sk, _, pk, _) = expand_decap_key::<L, MU>(&sk_seed);

            let mut enc_seed = [0u8; 32];
            let mut msg = [0u8; 32];
            rng.fill_bytes(&mut enc_seed);
            rng.fill_bytes(&mut msg);
            let ct_buf = &mut backing_buf[..T * RING_DEG / 8 + L * MODULUS_P_BITS * RING_DEG / 8];

            encrypt_deterministic::<L, MU, T>(&pk, &msg, &enc_seed, ct_buf);
            let recovered_msg = decrypt::<L, T>(&sk, ct_buf);
            assert_eq!(msg, recovered_msg);
        }
    }
}
