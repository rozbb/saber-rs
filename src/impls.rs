//! Contains external-facing impls of KEM traits

use crate::{
    consts::*,
    kem::{KemPublicKey, KemSecretKey},
    pke::ciphertext_len,
};

use core::convert::Infallible;

use kem_traits::{Decapsulate, Encapsulate};
use rand_core::CryptoRngCore;

/// A shared secret of a KEM execution
pub type SharedSecret = [u8; 32];

/// A secret key for the Saber KEM
pub struct SaberSecretKey(KemSecretKey<SABER_L>);

/// A public key for the Saber KEM
pub struct SaberPublicKey(KemPublicKey<SABER_L>);

/// A ciphertext, or "encapsulated key"", for the Saber KEM. This is a thin wrapper around
/// a fixed-size byte array.
#[derive(Clone)]
pub struct SaberCiphertext([u8; ciphertext_len::<SABER_L, SABER_MODULUS_T_BITS>()]);

impl SaberCiphertext {
    /// The length of the ciphertext, which is a byte array
    pub const LEN: usize = ciphertext_len::<SABER_L, SABER_MODULUS_T_BITS>();
}

impl SaberSecretKey {
    /// The length of the secret key when serialized to bytes
    pub const SERIALIZED_LEN: usize = KemSecretKey::<SABER_L>::SERIALIZED_LEN;

    /// Generate a fresh secret key
    pub fn generate(rng: &mut impl CryptoRngCore) -> Self {
        Self(KemSecretKey::generate::<SABER_MU>(rng))
    }

    /// Serialize this secret key into `out_buf`. `out_buf` MUST have length `Self::SERIALIZED_LEN`.
    ///
    /// # Panics
    /// Panics if `out_buf.len() != Self::SERIALIZED_LEN`.
    pub fn to_bytes(&self, out_buf: &mut [u8]) {
        self.0.to_bytes(out_buf);
    }

    /// Deserializes a secret key from `bytes`. `bytes` MUST have length `Self::SERIALIZED_LEN`.
    ///
    /// # Panics
    /// Panics if `bytes.len() != Self::SERIALIZED_LEN`.
    pub fn from_bytes(bytes: &[u8]) -> Self {
        Self(KemSecretKey::from_bytes(bytes))
    }

    /// Returns the public key corresponding to this secret key
    pub fn public_key(&self) -> SaberPublicKey {
        SaberPublicKey(self.0.public_key().clone())
    }
}

impl SaberPublicKey {
    /// The length of the public key when serialized to bytes
    pub const SERIALIZED_LEN: usize = KemPublicKey::<SABER_L>::SERIALIZED_LEN;

    /// Serialize this public key into `out_buf`. `out_buf` MUST have length `Self::SERIALIZED_LEN`.
    ///
    /// # Panics
    /// Panics if `out_buf.len() != Self::SERIALIZED_LEN`.
    pub fn to_bytes(&self, out_buf: &mut [u8]) {
        self.0.to_bytes(out_buf);
    }

    /// Deserializes a public key from `bytes`. `bytes` MUST have length `Self::SERIALIZED_LEN`.
    ///
    /// # Panics
    /// Panics if `bytes.len() != Self::SERIALIZED_LEN`.
    pub fn from_bytes(bytes: &[u8]) -> Self {
        Self(KemPublicKey::from_bytes(bytes))
    }
}

impl Default for SaberCiphertext {
    fn default() -> Self {
        SaberCiphertext([0u8; ciphertext_len::<SABER_L, SABER_MODULUS_T_BITS>()])
    }
}

impl AsRef<[u8]> for SaberCiphertext {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for SaberCiphertext {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl Encapsulate<SaberCiphertext, SharedSecret> for SaberPublicKey {
    /// Encapsulation cannot fail
    type Error = Infallible;

    fn encapsulate(
        &self,
        rng: &mut impl CryptoRngCore,
    ) -> Result<(SaberCiphertext, SharedSecret), Self::Error> {
        let mut ciphertext = SaberCiphertext::default();
        let shared_secret = crate::kem::encap::<SABER_L, SABER_MU, SABER_MODULUS_T_BITS>(
            rng,
            &self.0,
            &mut ciphertext.0,
        );

        Ok((ciphertext, shared_secret))
    }
}
impl Decapsulate<SaberCiphertext, SharedSecret> for SaberSecretKey {
    /// Decapsulation cannot fail
    type Error = Infallible;

    /// Decapsulates an encapsulated key and returns the resulting shared secret.
    /// If the encapsulated key is invalid, then the shared secret will be psuedorandom garbage.
    fn decapsulate(&self, encapsulated_key: &SaberCiphertext) -> Result<SharedSecret, Self::Error> {
        Ok(
            crate::kem::decap::<SABER_L, SABER_MU, SABER_MODULUS_T_BITS>(
                &self.0,
                &encapsulated_key.0,
            ),
        )
    }
}

#[test]
fn test_api() {
    let mut rng = rand::thread_rng();
    let sk = SaberSecretKey::generate(&mut rng);
    let pk = sk.public_key();

    // Serialize and deserialize the keys
    let mut sk_bytes = [0u8; SaberSecretKey::SERIALIZED_LEN];
    sk.to_bytes(&mut sk_bytes);
    let sk = SaberSecretKey::from_bytes(&sk_bytes);

    let mut pk_bytes = [0u8; SaberPublicKey::SERIALIZED_LEN];
    pk.to_bytes(&mut pk_bytes);
    let pk = SaberPublicKey::from_bytes(&pk_bytes);

    let (ct, ss1) = pk.encapsulate(&mut rng).unwrap();
    let ct_bytes = ct.as_ref();

    let mut receiver_ct = SaberCiphertext::default();
    receiver_ct.as_mut().copy_from_slice(ct_bytes);
    let ss2 = sk.decapsulate(&receiver_ct).unwrap();

    assert_eq!(ss1, ss2);
}