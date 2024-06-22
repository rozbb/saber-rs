//! Contains external-facing impls of KEM traits

use crate::{
    consts::*,
    kem::{KemPublicKey, KemSecretKey},
    pke::ciphertext_len,
};

use kem_traits::{Decapsulate, Encapsulate};
use rand_core::CryptoRngCore;

/// The shared secret of a KEM execution
pub type SharedSecret = [u8; 32];

/// The private key for the Saber KEM
pub struct SaberPrivkey(KemSecretKey<SABER_L>);

/// The public key for the Saber KEM
pub struct SaberPubkey(KemPublicKey<SABER_L>);

/// The ciphertext, or "encapsulated key"", for the Saber KEM
#[derive(Clone)]
pub struct SaberCiphertext([u8; ciphertext_len::<SABER_L, SABER_MODULUS_T_BITS>()]);

impl SaberPrivkey {
    pub const SERIALIZED_LEN: usize = KemSecretKey::<SABER_L>::SERIALIZED_LEN;

    /// Generate a fresh private key
    pub fn generate(rng: &mut impl CryptoRngCore) -> Self {
        Self(KemSecretKey::generate::<SABER_MU>(rng))
    }

    pub fn to_bytes(&self, out_buf: &mut [u8]) {
        self.0.to_bytes(out_buf);
    }

    pub fn from_bytes(bytes: &[u8]) -> Self {
        Self(KemSecretKey::from_bytes(bytes))
    }

    pub fn public_key(&self) -> SaberPubkey {
        SaberPubkey(self.0.public_key().clone())
    }
}

impl SaberPubkey {
    pub const SERIALIZED_LEN: usize = KemPublicKey::<SABER_L>::SERIALIZED_LEN;

    pub fn to_bytes(&self, out_buf: &mut [u8]) {
        self.0.to_bytes(out_buf);
    }

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

impl Encapsulate<SaberCiphertext, SharedSecret> for SaberPubkey {
    type Error = ();

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
impl Decapsulate<SaberCiphertext, SharedSecret> for SaberPrivkey {
    type Error = ();

    fn decapsulate(&self, encapsulated_key: &SaberCiphertext) -> Result<SharedSecret, ()> {
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
    let sk = SaberPrivkey::generate(&mut rng);
    let pk = sk.public_key();

    // Serialize and deserialize the keys
    let mut sk_bytes = [0u8; SaberPrivkey::SERIALIZED_LEN];
    sk.to_bytes(&mut sk_bytes);
    let sk = SaberPrivkey::from_bytes(&sk_bytes);

    let mut pk_bytes = [0u8; SaberPubkey::SERIALIZED_LEN];
    pk.to_bytes(&mut pk_bytes);
    let pk = SaberPubkey::from_bytes(&pk_bytes);

    let (ct, ss1) = pk.encapsulate(&mut rng).unwrap();
    let ct_bytes = ct.as_ref();

    let mut receiver_ct = SaberCiphertext::default();
    receiver_ct.as_mut().copy_from_slice(ct_bytes);
    let ss2 = sk.decapsulate(&receiver_ct).unwrap();

    assert_eq!(ss1, ss2);
}
