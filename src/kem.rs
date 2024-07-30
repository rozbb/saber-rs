//! This file implements the IND-CCA-secure Saber.KEM scheme

use crate::pke::{
    self, ciphertext_len, max_ciphertext_len, max_pke_pubkey_serialized_len, PkePublicKey,
    PkeSecretKey,
};

use rand_core::CryptoRngCore;
use sha3::{digest::Digest, Sha3_256};
use subtle::{ConditionallySelectable, ConstantTimeEq};

/// A public key for the IND-CCA-secure Saber KEM scheme
pub type KemPublicKey<const L: usize> = PkePublicKey<L>;

/// The shared secret of a KEM operation
pub type SharedSecret = [u8; 32];

/// A secret key for the IND-CCA-secure Saber KEM scheme
pub struct KemSecretKey<const L: usize> {
    /// Used for deriving pseudorandom shared secrets when decap fails
    z: [u8; 32],
    /// The PKE secret key
    pke_sk: PkeSecretKey<L>,
    /// The PKE public key that `cpa_sk` generated
    pke_pk: PkePublicKey<L>,
    /// The hash of `cpa_pk`
    hash_pke_pk: [u8; 32],
}

impl<const L: usize> KemSecretKey<L> {
    pub const SERIALIZED_LEN: usize =
        PkeSecretKey::<L>::SERIALIZED_LEN + PkePublicKey::<L>::SERIALIZED_LEN + 32 + 32;

    // Algorithm 20, Saber.KEM.KeyGen
    /// Generate a fresh secret key
    pub fn generate<const MU: usize>(rng: &mut impl CryptoRngCore) -> KemSecretKey<L> {
        let (pke_sk, pke_pk) = pke::gen_keypair::<L, MU>(rng);

        // Hash the public key
        let hash_pke_pk = {
            // Serialize the pubkey into a buffer of the appropriate size
            let mut buf = [0u8; max_pke_pubkey_serialized_len()];
            let pk_bytes = &mut buf[..PkePublicKey::<L>::SERIALIZED_LEN];
            pke_pk.to_bytes(pk_bytes);
            Sha3_256::digest(pk_bytes)
        };

        // Sample some random bytes for z
        let mut z = [0u8; 32];
        rng.fill_bytes(&mut z);

        KemSecretKey {
            z,
            hash_pke_pk: hash_pke_pk.into(),
            pke_pk,
            pke_sk,
        }
    }

    /// Serialize this secret key into `out_buf`. `out_buf` MUST have length `Self::SERIALIZED_LEN`
    pub fn to_bytes(&self, out_buf: &mut [u8]) {
        assert_eq!(out_buf.len(), Self::SERIALIZED_LEN);
        let rest = out_buf;

        // Serialization order, according to reference impl, is sk, pk, hash_pk, z
        // https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/kem.c#L12-L25

        let (out_cpa_sk, rest) = rest.split_at_mut(PkeSecretKey::<L>::SERIALIZED_LEN);
        self.pke_sk.to_bytes(out_cpa_sk);

        let (out_cpa_pk, rest) = rest.split_at_mut(PkePublicKey::<L>::SERIALIZED_LEN);
        self.pke_pk.to_bytes(out_cpa_pk);

        let (out_hash_pk, rest) = rest.split_at_mut(32);
        out_hash_pk.copy_from_slice(&self.hash_pke_pk);

        let (out_z, rest) = rest.split_at_mut(32);
        out_z.copy_from_slice(&self.z);

        assert_eq!(rest.len(), 0);
    }

    /// Deserializes a `KemSecretKey` from the given byte string. `bytes` MUST have length
    /// `Self::SERIALIZED_LEN`.
    pub fn from_bytes(bytes: &[u8]) -> Self {
        assert_eq!(bytes.len(), Self::SERIALIZED_LEN);

        // Serialization order is sk, pk, hash_pk, z

        let (bytes, rest) = bytes.split_at(PkeSecretKey::<L>::SERIALIZED_LEN);
        let pke_sk = PkeSecretKey::from_bytes(bytes);

        let (bytes, rest) = rest.split_at(PkePublicKey::<L>::SERIALIZED_LEN);
        let pke_pk = PkePublicKey::from_bytes(bytes);

        let (bytes, rest) = rest.split_at(32);
        let hash_pke_pk = bytes.try_into().expect("split_at 32 is 32-bytes long");

        let (z, rest) = rest.split_at(32);
        let z = z.try_into().expect("split_at 32 is 32-bytes long");

        assert_eq!(rest.len(), 0);

        Self {
            z,
            hash_pke_pk,
            pke_pk,
            pke_sk,
        }
    }

    pub(crate) fn public_key(&self) -> &KemPublicKey<L> {
        &self.pke_pk
    }
}

// Algorithm 21, Saber.KEM.Encaps
/// Encapsulate a shared secret to the given public key. Returns the shared secret.
/// `out_buf` MUST have length `ciphertext_len::<L>()`.
pub(crate) fn encap<const L: usize, const MU: usize, const MODULUS_T_BITS: usize>(
    rng: &mut impl CryptoRngCore,
    kem_pk: &KemPublicKey<L>,
    out_buf: &mut [u8],
) -> SharedSecret {
    // Before using the message, we need to hash it
    let mut m_unhashed = [0u8; 32];
    rng.fill_bytes(&mut m_unhashed);
    let m = Sha3_256::digest(m_unhashed).into();

    encap_deterministic::<L, MU, MODULUS_T_BITS>(&m, kem_pk, out_buf)
}

// Algorithm 21, Saber.KEM.Encaps (cont.)
/// Encapsulate a shared secret to the given public key. Uses the given message `m` as the data
/// being encrypted to `kem_pk` Returns the shared secret. `out_buf` MUST have length
/// `ciphertext_len::<L>()`.
pub(crate) fn encap_deterministic<const L: usize, const MU: usize, const MODULUS_T_BITS: usize>(
    m: &[u8; 32],
    kem_pk: &KemPublicKey<L>,
    out_buf: &mut [u8],
) -> SharedSecret {
    // The PKE public key is the same as the KEM public key
    let pke_pk = kem_pk;

    // Hash the public key
    let hash_pke_pk = {
        // Serialize the pubkey into a buffer of the appropriate size
        let mut buf = [0u8; max_pke_pubkey_serialized_len()];
        let pk_bytes = &mut buf[..PkePublicKey::<L>::SERIALIZED_LEN];
        pke_pk.to_bytes(pk_bytes);
        Sha3_256::digest(pk_bytes)
    };

    // kr = SHA3-512(m || hash_pke_pk)
    // The spec has the hash input order switched, but we're following the reference impl
    // https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/kem.c#L35-L39
    let kr = sha3::Sha3_512::new()
        .chain_update(m)
        .chain_update(hash_pke_pk)
        .finalize();

    // k || r = kr
    // The spec has the order switched, but, again, we're following the reference impl
    // https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/kem.c#L42-L46
    let (k, r) = kr.split_at(32);
    let k: [u8; 32] = k.try_into().expect("split-at 32 is 32-bytes long");
    let r: [u8; 32] = r
        .try_into()
        .expect("[u8; 64].split_at(32).1 is 32-bytes long");

    // ct = Saber.PKE.Enc_pk(m; r)
    pke::encrypt_deterministic::<L, MU, MODULUS_T_BITS>(kem_pk, m, &r, out_buf);

    // r' = SHA3-256(ct)
    let rprime = Sha3_256::digest(out_buf);

    // shared secret = SHA3-256(k || r')
    // The spec has the hash input order switched, but we're following the reference impl
    // https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/kem.c#L46

    Sha3_256::new()
        .chain_update(k)
        .chain_update(rprime)
        .finalize()
        .into()
}

// Algorithm 22, Saber.KEM.Decaps
/// Decapsulates a shared secret from the given ciphertext and secret key.
/// Returns the shared secret or a random value if the ciphertext is invalid.
/// `ciphertext` MUST have length `ciphertext_len::<L>()`.
pub fn decap<const L: usize, const MU: usize, const MODULUS_T_BITS: usize>(
    sk: &KemSecretKey<L>,
    ciphertext: &[u8],
) -> SharedSecret {
    assert_eq!(ciphertext.len(), ciphertext_len::<L, MODULUS_T_BITS>());

    // Decrypt the PKE ciphertext
    let m = pke::decrypt::<L, MODULUS_T_BITS>(&sk.pke_sk, ciphertext);

    // kr = SHA3-512(m || hash_pk)
    // The spec has the hash input order switched, but we're following the reference impl
    // https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/kem.c#L35-L39
    let kr = sha3::Sha3_512::new()
        .chain_update(m)
        .chain_update(sk.hash_pke_pk)
        .finalize();

    // k || r = kr
    // The spec has the order switched, but, again, we're following the reference impl
    // https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/kem.c#L42-L46
    let (k, r) = kr.split_at(32);
    let k: [u8; 32] = k.try_into().expect("split-at 32 is 32-bytes long");
    let r: [u8; 32] = r
        .try_into()
        .expect("[u8; 64].split_at(32).1 is 32-bytes long");

    // ct' = Saber.PKE.Enc_pk(m; r)
    let mut buf = [0u8; max_ciphertext_len()];
    let reconstructed_ct = &mut buf[..ciphertext_len::<L, MODULUS_T_BITS>()];
    pke::encrypt_deterministic::<L, MU, MODULUS_T_BITS>(&sk.pke_pk, &m, &r, reconstructed_ct);

    // r' = SHA3-256(ct)
    let rprime = Sha3_256::digest(ciphertext);

    // k_or_z = k if reconstruction matched, else z. We do this in constant time using `subtle`
    let reconstruction_matched = reconstructed_ct.ct_eq(ciphertext);
    let k_or_z = <[u8; 32]>::conditional_select(&sk.z, &k, reconstruction_matched);

    // session key = SHA3-256(k_or_z || r')
    // The spec has the hash input order switched, but we're following the reference impl
    // https://github.com/KULeuven-COSIC/SABER/blob/f7f39e4db2f3e22a21e1dd635e0601caae2b4510/Reference_Implementation_KEM/kem.c#L46

    Sha3_256::new()
        .chain_update(k_or_z)
        .chain_update(rprime)
        .finalize()
        .into()
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::consts::*;

    use rand::Rng;

    #[test]
    fn lightsaber_cca_kem() {
        test_encap_decap::<LIGHTSABER_L, LIGHTSABER_MU, LIGHTSABER_MODULUS_T_BITS>();
    }

    #[test]
    fn saber_cca_kem() {
        test_encap_decap::<SABER_L, SABER_MU, SABER_MODULUS_T_BITS>();
    }

    #[test]
    fn firesaber_cca_kem() {
        test_encap_decap::<FIRESABER_L, FIRESABER_MU, FIRESABER_MODULUS_T_BITS>();
    }

    fn test_encap_decap<const L: usize, const MU: usize, const MODULUS_T_BITS: usize>() {
        let mut rng = rand::thread_rng();
        let mut backing_buf = [0u8; max_ciphertext_len()];

        for _ in 0..100 {
            let sk = KemSecretKey::<L>::generate::<MU>(&mut rng);
            let pk = sk.public_key();
            let ct_buf = &mut backing_buf[..ciphertext_len::<L, MODULUS_T_BITS>()];

            let m: [u8; 32] = rng.gen();
            let ss1 = encap_deterministic::<L, MU, MODULUS_T_BITS>(&m, pk, ct_buf);
            let ss2 = decap::<L, MU, MODULUS_T_BITS>(&sk, &ct_buf);
            assert_eq!(ss1, ss2);

            // Check that the Fujisaki-Okamoto transform was implemented properly. That is, a
            // perturbed ciphertext should yield a secret key that is totally unguessable to the
            // encapsulator
            let perturbed_ct = ct_buf;
            // XOR the ciphertext with a random (nonzero) byte in a random location
            let idx = (rng.gen::<u32>() as usize) % perturbed_ct.len();
            let byte = loop {
                let b = rng.gen::<u8>();
                if b != 0 {
                    break b;
                }
            };
            perturbed_ct[idx] ^= byte;
            // Decapsulate the perturbed ciphertext
            let ss1 = decap::<L, MU, MODULUS_T_BITS>(&sk, &perturbed_ct);
            // Try to guess what the decapsulation would be. This is the value that it would be if
            // we messed up the F-O transform and accidentally returned SHA3-256(k || r')
            // regardless of whether the equality check succeeded
            let ss2 = recompute_shared_secret_for_pertrubed_ciphertext::<L, MU, MODULUS_T_BITS>(
                &m,
                pk,
                &perturbed_ct,
            );
            assert_ne!(ss1, ss2);
        }
    }

    /// Model an adversary who has encapsulated a value to a given public key and has perturbed the
    /// ciphertext a little. If Fujisaki-Okamoto is not implemented, they are able to recalculate
    /// the new shared secret without knowledge of the secret key
    pub(crate) fn recompute_shared_secret_for_pertrubed_ciphertext<
        const L: usize,
        const MU: usize,
        const MODULUS_T_BITS: usize,
    >(
        m: &[u8; 32],
        kem_pk: &KemPublicKey<L>,
        ct: &[u8],
    ) -> SharedSecret {
        // Do all calculations as before
        let pke_pk = kem_pk;
        let hash_pke_pk = {
            // Serialize the pubkey into a buffer of the appropriate size
            let mut buf = [0u8; max_pke_pubkey_serialized_len()];
            let pk_bytes = &mut buf[..PkePublicKey::<L>::SERIALIZED_LEN];
            pke_pk.to_bytes(pk_bytes);
            Sha3_256::digest(pk_bytes)
        };
        let kr = sha3::Sha3_512::new()
            .chain_update(m)
            .chain_update(hash_pke_pk)
            .finalize();
        let (k, _) = kr.split_at(32);

        // Only difference: use the given ciphertext to derive r'
        let rprime = Sha3_256::digest(ct);

        // Return SHA3-256(k || r')
        Sha3_256::new()
            .chain_update(k)
            .chain_update(rprime)
            .finalize()
            .into()
    }
}
