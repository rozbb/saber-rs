//! This file implements the IND-CCA-secure Kopis KEM scheme

use crate::{
    consts::{DOMSEP_FO, DOMSEP_NOREJECT, DOMSEP_PKHASH},
    pke::{self, ciphertext_len, expand_decap_key, max_ciphertext_len, PkePublicKey, PkeSecretKey},
    turboshake256_hash,
};

use rand_core::CryptoRng;
use subtle::{ConditionallySelectable, ConstantTimeEq};
use turboshake::digest::{ExtendableOutput, Update, XofReader};
use turboshake::CTurboShake256;

/// A public key for the IND-CCA-secure Kopis KEM scheme
pub struct KemPublicKey<const L: usize> {
    /// The PKE public key
    pke_pk: PkePublicKey<L>,
    /// The hash of `pke_pk`
    hash_pke_pk: [u8; 32],
}

impl<const L: usize> KemPublicKey<L> {
    pub(crate) const SERIALIZED_LEN: usize = PkePublicKey::<L>::SERIALIZED_LEN;

    /// Serializes just `pke_pk`
    pub(crate) fn to_bytes(&self, out_buf: &mut [u8]) {
        self.pke_pk.to_bytes(out_buf);
    }

    /// Deserializes from `pke_pk`, and recomputes `hash_pke_pk`
    pub(crate) fn from_bytes(bytes: &[u8]) -> Self {
        let pke_pk = PkePublicKey::from_bytes(bytes);
        // Recompute the hash: pkh = TurboSHAKE256(pk, 32, DOMSEP_PKHASH)
        let hash_pke_pk = turboshake256_hash::<DOMSEP_PKHASH>(&[bytes]);

        KemPublicKey {
            pke_pk,
            hash_pke_pk,
        }
    }
}

/// The shared secret of a KEM operation
pub type SharedSecret = [u8; 32];

/// A secret key for the IND-CCA-secure Kopis KEM scheme.
///
/// The canonical secret key is a 32-byte seed. This struct stores the expanded form for
/// efficiency (avoiding re-expansion on every decapsulation).
pub struct KemSecretKey<const L: usize> {
    /// The 32-byte seed (the canonical secret key, used for serialization)
    seed: [u8; 32],
    /// Used for deriving pseudorandom shared secrets when decap fails
    z: [u8; 32],
    /// The PKE secret key (expanded from seed)
    pke_sk: PkeSecretKey<L>,
    /// The PKE public key that `pke_sk` generated
    pke_pk: PkePublicKey<L>,
    /// The hash of `pke_pk`
    hash_pke_pk: [u8; 32],
}

impl<const L: usize> KemSecretKey<L> {
    /// Construct a secret key from a 32-byte seed by expanding it via `ExpandDecapKey`.
    pub fn from_bytes<const MU: usize>(seed: &[u8; 32]) -> KemSecretKey<L> {
        let (pke_sk, z, pke_pk, hash_pke_pk) = expand_decap_key::<L, MU>(seed);

        KemSecretKey {
            seed: *seed,
            z,
            pke_sk,
            pke_pk,
            hash_pke_pk,
        }
    }

    /// Generate a fresh secret key by sampling a random 32-byte seed.
    pub fn generate<const MU: usize>(rng: &mut impl CryptoRng) -> KemSecretKey<L> {
        let mut seed = [0u8; 32];
        rng.fill_bytes(&mut seed);
        Self::from_bytes::<MU>(&seed)
    }

    /// Serialize this secret key to bytes. The secret key is just the 32-byte seed.
    pub fn to_bytes(&self) -> [u8; 32] {
        self.seed
    }

    pub(crate) fn public_key(&self) -> KemPublicKey<L> {
        KemPublicKey {
            pke_pk: self.pke_pk.clone(),
            hash_pke_pk: self.hash_pke_pk,
        }
    }
}

/// Encapsulate a shared secret to the given public key. Returns the shared secret.
/// `out_buf` MUST have length `ciphertext_len::<L, T>()`.
pub(crate) fn encap<const L: usize, const MU: usize, const T: usize>(
    rng: &mut impl CryptoRng,
    kem_pk: &KemPublicKey<L>,
    out_buf: &mut [u8],
) -> SharedSecret {
    let mut randomness = [0u8; 32];
    rng.fill_bytes(&mut randomness);

    encap_deterministic::<L, MU, T>(&randomness, kem_pk, out_buf)
}

/// Encapsulate a shared secret to the given public key using the given `randomness`.
/// Returns the shared secret. `out_buf` MUST have length `ciphertext_len::<L, T>()`.
pub(crate) fn encap_deterministic<const L: usize, const MU: usize, const T: usize>(
    randomness: &[u8; 32],
    kem_pk: &KemPublicKey<L>,
    out_buf: &mut [u8],
) -> SharedSecret {
    let KemPublicKey {
        pke_pk,
        hash_pke_pk,
    } = kem_pk;

    // k || r = TurboSHAKE256()
    let mut k = [0u8; 32];
    let mut r = [0u8; 32];

    let mut xof = {
        let mut hasher = CTurboShake256::<DOMSEP_FO>::default();
        hasher.update(randomness);
        hasher.update(hash_pke_pk);
        hasher.finalize_xof()
    };
    xof.read(&mut k);
    xof.read(&mut r);

    // c = PkeEncrypt(r, pk, randomness)
    pke::encrypt_deterministic::<L, MU, T>(pke_pk, randomness, &r, out_buf);

    k
}

/// Decapsulates a shared secret from the given ciphertext and secret key.
/// Returns the shared secret or a pseudorandom value if the ciphertext is invalid.
/// `ciphertext` MUST have length `ciphertext_len::<L, T>()`.
pub fn decap<const L: usize, const MU: usize, const T: usize>(
    sk: &KemSecretKey<L>,
    ciphertext: &[u8],
) -> SharedSecret {
    assert_eq!(ciphertext.len(), ciphertext_len::<L, T>());

    // randomness = PkeDecrypt(sk, c)
    let randomness = pke::decrypt::<L, T>(&sk.pke_sk, ciphertext);

    // k || rprime = TurboSHAKE256()
    let mut k = [0u8; 32];
    let mut rprime = [0u8; 32];

    let mut xof = {
        let mut hasher = CTurboShake256::<DOMSEP_FO>::default();
        hasher.update(&randomness);
        hasher.update(&sk.hash_pke_pk);
        hasher.finalize_xof()
    };
    xof.read(&mut k);
    xof.read(&mut rprime);

    // cprime = PkeEncrypt(rprime, pk, randomness)
    let mut buf = [0u8; max_ciphertext_len()];
    let reconstructed_ct = &mut buf[..ciphertext_len::<L, T>()];
    pke::encrypt_deterministic::<L, MU, T>(&sk.pke_pk, &randomness, &rprime, reconstructed_ct);

    // Compute rejection value: TurboSHAKE256(z || c, 32, DOMSEP_NOREJECT)
    let reject_val = turboshake256_hash::<DOMSEP_NOREJECT>(&[&sk.z, ciphertext]);

    // Constant-time select: return k if c == cprime, else return reject_val
    let reconstruction_matched = reconstructed_ct.ct_eq(ciphertext);
    <[u8; 32]>::conditional_select(&reject_val, &k, reconstruction_matched)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::consts::*;
    use crate::pke::max_pke_pubkey_serialized_len;

    use rand::Rng;

    #[test]
    fn kopis512_cca_kem() {
        test_encap_decap::<KOPIS512_L, KOPIS512_MU, KOPIS512_T>();
    }

    #[test]
    fn kopis768_cca_kem() {
        test_encap_decap::<KOPIS768_L, KOPIS768_MU, KOPIS768_T>();
    }

    #[test]
    fn kopis1024_cca_kem() {
        test_encap_decap::<KOPIS1024_L, KOPIS1024_MU, KOPIS1024_T>();
    }

    fn test_encap_decap<const L: usize, const MU: usize, const T: usize>() {
        let mut rng = rand::rng();
        let mut backing_buf = [0u8; max_ciphertext_len()];

        for _ in 0..100 {
            let sk = KemSecretKey::<L>::generate::<MU>(&mut rng);
            let pk = sk.public_key();
            let ct_buf = &mut backing_buf[..ciphertext_len::<L, T>()];

            let randomness: [u8; 32] = rng.random();
            let ss1 = encap_deterministic::<L, MU, T>(&randomness, &pk, ct_buf);
            let ss2 = decap::<L, MU, T>(&sk, ct_buf);
            assert_eq!(ss1, ss2);

            // Check that the Fujisaki-Okamoto transform was implemented properly. That is, a
            // perturbed ciphertext should yield a secret key that is totally unguessable to the
            // encapsulator
            let perturbed_ct = ct_buf;
            // XOR the ciphertext with a random (nonzero) byte in a random location
            let idx = (rng.random::<u32>() as usize) % perturbed_ct.len();
            let byte = loop {
                let b = rng.random::<u8>();
                if b != 0 {
                    break b;
                }
            };
            perturbed_ct[idx] ^= byte;
            // Decapsulate the perturbed ciphertext
            let ss1 = decap::<L, MU, T>(&sk, perturbed_ct);
            // Try to guess what the decapsulation would be. If we messed up the F-O transform
            // and returned k regardless of the equality check, the adversary (who knows
            // randomness and pk) could compute k.
            let ss2 =
                recompute_shared_secret_for_pertrubed_ciphertext::<L, MU, T>(&randomness, &pk);
            assert_ne!(ss1, ss2);
        }
    }

    /// Model an adversary who has encapsulated a value to a given public key and has perturbed the
    /// ciphertext. If Fujisaki-Okamoto is not implemented, they could predict the shared secret
    /// by recomputing k from the randomness and public key hash.
    pub(crate) fn recompute_shared_secret_for_pertrubed_ciphertext<
        const L: usize,
        const MU: usize,
        const T: usize,
    >(
        randomness: &[u8; 32],
        kem_pk: &KemPublicKey<L>,
    ) -> SharedSecret {
        // Recompute pkh
        let hash_pke_pk = {
            let mut buf = [0u8; max_pke_pubkey_serialized_len()];
            let pk_bytes = &mut buf[..PkePublicKey::<L>::SERIALIZED_LEN];
            kem_pk.to_bytes(pk_bytes);

            turboshake256_hash::<DOMSEP_PKHASH>(&[pk_bytes])
        };

        // k || r = TurboSHAKE256(randomness || pkh, 64, DOMSEP_FO)
        let mut k = [0u8; 32];

        let mut xof = {
            let mut hasher = CTurboShake256::<DOMSEP_FO>::default();
            hasher.update(randomness);
            hasher.update(&hash_pke_pk);
            hasher.finalize_xof()
        };
        xof.read(&mut k);

        k
    }
}
