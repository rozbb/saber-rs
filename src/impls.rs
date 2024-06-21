//! Contains external-facing impls of KEM traits

use crate::kem::KemSecretKey;
use crate::pke::PkePublicKey;
use crate::{consts::*, pke::ciphertext_len};

use kem_traits::{Decapsulate, Encapsulate};

use rand_core::CryptoRngCore;

pub type SaberPrivkey = KemSecretKey<SABER_L>;
pub type SaberPubkey = PkePublicKey<SABER_L>;
pub type SaberCiphertext = [u8; ciphertext_len::<SABER_L, SABER_MODULUS_T_BITS>()];
pub type SaberSharedSecret = [u8; 32];

impl Encapsulate<SaberCiphertext, SaberSharedSecret> for SaberPubkey {
    type Error = ();

    fn encapsulate(
        &self,
        rng: &mut impl CryptoRngCore,
    ) -> Result<(SaberCiphertext, SaberSharedSecret), Self::Error> {
        let mut ciphertext = [0u8; ciphertext_len::<SABER_L, SABER_MODULUS_T_BITS>()];
        let shared_secret = crate::kem::encap::<SABER_L, SABER_MU, SABER_MODULUS_T_BITS>(
            rng,
            self,
            &mut ciphertext,
        );

        Ok((ciphertext, shared_secret))
    }
}
impl Decapsulate<SaberCiphertext, SaberSharedSecret> for SaberPrivkey {
    type Error = ();

    fn decapsulate(&self, encapsulated_key: &SaberCiphertext) -> Result<SaberSharedSecret, ()> {
        Ok(crate::kem::decap::<SABER_L, SABER_MU, SABER_MODULUS_T_BITS>(&self, encapsulated_key))
    }
}
