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

/// Defines convenience types and impls for a given Saber variant
macro_rules! variant_impl {
    (
        $variant_name:ident,
        $mod_doc:expr,
        $pubkey_name:ident,
        $privkey_name:ident,
        $ciphertext_name:ident,
        $variant_ell:expr,
        $variant_mu:expr,
        $variant_modt_bits:expr
    ) => {
        #[doc = $mod_doc]
        pub mod $variant_name {
            use super::*;

            /// A secret key for this KEM
            pub struct $privkey_name(KemSecretKey<$variant_ell>);

            /// A public key for this KEM
            pub struct $pubkey_name(KemPublicKey<$variant_ell>);

            /// A ciphertext, or "encapsulated key"", for this KEM. This is a thin wrapper around
            /// a fixed-size byte array.
            #[derive(Clone)]
            pub struct $ciphertext_name([u8; ciphertext_len::<$variant_ell, $variant_modt_bits>()]);

            impl $ciphertext_name {
                /// The length of the ciphertext, which is a byte array
                pub const LEN: usize = ciphertext_len::<$variant_ell, $variant_modt_bits>();
            }

            impl $privkey_name {
                /// The length of the secret key when serialized to bytes
                pub const SERIALIZED_LEN: usize = KemSecretKey::<$variant_ell>::SERIALIZED_LEN;

                /// Generate a fresh secret key
                pub fn generate(rng: &mut impl CryptoRngCore) -> Self {
                    Self(KemSecretKey::generate::<$variant_mu>(rng))
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
                pub fn public_key(&self) -> $pubkey_name {
                    $pubkey_name(self.0.public_key().clone())
                }
            }

            impl $pubkey_name {
                /// The length of the public key when serialized to bytes
                pub const SERIALIZED_LEN: usize = KemPublicKey::<$variant_ell>::SERIALIZED_LEN;

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

            impl Default for $ciphertext_name {
                fn default() -> Self {
                    $ciphertext_name([0u8; ciphertext_len::<$variant_ell, $variant_modt_bits>()])
                }
            }

            impl AsRef<[u8]> for $ciphertext_name {
                fn as_ref(&self) -> &[u8] {
                    &self.0
                }
            }

            impl AsMut<[u8]> for $ciphertext_name {
                fn as_mut(&mut self) -> &mut [u8] {
                    &mut self.0
                }
            }

            impl Encapsulate<$ciphertext_name, SharedSecret> for $pubkey_name {
                /// Encapsulation cannot fail
                type Error = Infallible;

                /// Encapsulates a fresh shared secret. This cannot fail.
                fn encapsulate(
                    &self,
                    rng: &mut impl CryptoRngCore,
                ) -> Result<($ciphertext_name, SharedSecret), Infallible> {
                    let mut ciphertext = $ciphertext_name::default();
                    let shared_secret = crate::kem::encap::<
                        $variant_ell,
                        $variant_mu,
                        $variant_modt_bits,
                    >(rng, &self.0, &mut ciphertext.0);

                    Ok((ciphertext, shared_secret))
                }
            }

            impl Decapsulate<$ciphertext_name, SharedSecret> for $privkey_name {
                /// Decapsulation cannot fail
                type Error = Infallible;

                /// Decapsulates an encapsulated key and returns the resulting shared secret.
                /// If the encapsulated key is invalid, then the shared secret will be psuedorandom garbage.
                fn decapsulate(
                    &self,
                    encapsulated_key: &$ciphertext_name,
                ) -> Result<SharedSecret, Infallible> {
                    Ok(crate::kem::decap::<
                        $variant_ell,
                        $variant_mu,
                        $variant_modt_bits,
                    >(&self.0, &encapsulated_key.0))
                }
            }

            /// Basic test that keygen, encap, decap, ser, and deser work
            #[test]
            fn test_api() {
                let mut rng = rand::thread_rng();
                let sk = $privkey_name::generate(&mut rng);
                let pk = sk.public_key();

                // Serialize and deserialize the keys
                let mut sk_bytes = [0u8; $privkey_name::SERIALIZED_LEN];
                sk.to_bytes(&mut sk_bytes);
                let sk = $privkey_name::from_bytes(&sk_bytes);

                let mut pk_bytes = [0u8; $pubkey_name::SERIALIZED_LEN];
                pk.to_bytes(&mut pk_bytes);
                let pk = $pubkey_name::from_bytes(&pk_bytes);

                let (ct, ss1) = pk.encapsulate(&mut rng).unwrap();
                let ct_bytes = ct.as_ref();

                let mut receiver_ct = $ciphertext_name::default();
                receiver_ct.as_mut().copy_from_slice(ct_bytes);
                let ss2 = sk.decapsulate(&receiver_ct).unwrap();

                assert_eq!(ss1, ss2);
            }
        }
    };
}

variant_impl!(
    lightsaber,
    "LightSaber is designed to have security close to that of AES-128",
    LightsaberPublicKey,
    LightsaberSecretKey,
    LightsaberCiphertext,
    LIGHTSABER_L,
    LIGHTSABER_MU,
    LIGHTSABER_MODULUS_T_BITS
);

variant_impl!(
    saber,
    "Saber is designed to have security close to that of AES-192",
    SaberPublicKey,
    SaberSecretKey,
    SaberCiphertext,
    SABER_L,
    SABER_MU,
    SABER_MODULUS_T_BITS
);

variant_impl!(
    firesaber,
    "FireSaber is designed to have security close to that of AES-256",
    FiresaberPublicKey,
    FiresaberSecretKey,
    FiresaberCiphertext,
    FIRESABER_L,
    FIRESABER_MU,
    FIRESABER_MODULUS_T_BITS
);
