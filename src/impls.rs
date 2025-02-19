//! This file implements the external-facing API of our KEM

use crate::{
    consts::*,
    kem::{KemPublicKey, KemSecretKey},
    pke::ciphertext_len,
};

use rand_core::CryptoRng;
use zeroize::{Zeroize, ZeroizeOnDrop};

/// A shared secret of a KEM execution. This is just a `[u8; 32]` that zeroes itself from memory
/// when it goes out of scope.
#[derive(Zeroize, ZeroizeOnDrop)]
pub struct SharedSecret([u8; 32]);

impl SharedSecret {
    /// Returns the shared secret as a slice
    #[inline]
    pub fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

/// Defines convenience types and impls for a given Saber variant
macro_rules! variant_impl {
    (
        $variant_name:ident,
        $mod_doc:expr,
        $pubkey_name:ident,
        $privkey_name:ident,
        $ciphertext_name:ident,
        $ciphertext_len_name:ident,
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

            /// The length of a ciphertext, or "encapsulated key", for this KEM
            pub const $ciphertext_len_name: usize =
                ciphertext_len::<$variant_ell, $variant_modt_bits>();

            /// A ciphertext, or "encapsulated key", for this KEM. This is just a bytestring with
            /// length `
            #[doc = stringify!($ciphertext_len_name)]
            /// `.
            pub type $ciphertext_name = [u8; $ciphertext_len_name];

            impl $privkey_name {
                /// The length of the secret key when serialized to bytes
                pub const SERIALIZED_LEN: usize = KemSecretKey::<$variant_ell>::SERIALIZED_LEN;

                /// Generate a fresh secret key
                pub fn generate(rng: &mut impl CryptoRng) -> Self {
                    Self(KemSecretKey::generate::<$variant_mu>(rng))
                }

                /// Serializes this secret key into `out_buf`, of length `Self::SERIALIZED_LEN`
                pub fn to_bytes(&self, out_buf: &mut [u8; Self::SERIALIZED_LEN]) {
                    self.0.to_bytes(out_buf);
                }

                /// Deserializes a secret key from `bytes`, of length `Self::SERIALIZED_LEN`
                pub fn from_bytes(bytes: &[u8; Self::SERIALIZED_LEN]) -> Self {
                    Self(KemSecretKey::from_bytes(bytes))
                }

                /// Returns the public key corresponding to this secret key
                pub fn public_key(&self) -> $pubkey_name {
                    $pubkey_name(self.0.public_key())
                }
            }

            impl $pubkey_name {
                /// The length of the public key when serialized to bytes
                pub const SERIALIZED_LEN: usize = KemPublicKey::<$variant_ell>::SERIALIZED_LEN;

                /// Serializes this public key into `out_buf`, of length `Self::SERIALIZED_LEN`
                pub fn to_bytes(&self, out_buf: &mut [u8; Self::SERIALIZED_LEN]) {
                    self.0.to_bytes(out_buf);
                }

                /// Deserializes a public key from `bytes`, of length `Self::SERIALIZED_LEN`
                pub fn from_bytes(bytes: &[u8; Self::SERIALIZED_LEN]) -> Self {
                    Self(KemPublicKey::from_bytes(bytes))
                }
            }

            impl $pubkey_name {
                /// Encapsulates a fresh shared secret
                pub fn encapsulate(
                    &self,
                    rng: &mut impl CryptoRng,
                ) -> ($ciphertext_name, SharedSecret) {
                    let mut ct = [0u8; $ciphertext_len_name];
                    let ss = self.encapsulate_in_place(rng, &mut ct);
                    (ct, ss)
                }
            }

            impl $pubkey_name {
                /// Encapsulates a fresh shared secret and place the ciphertext in the given
                /// buffer
                pub fn encapsulate_in_place(
                    &self,
                    rng: &mut impl CryptoRng,
                    ct_out: &mut $ciphertext_name,
                ) -> SharedSecret {
                    let shared_secret =
                        crate::kem::encap::<$variant_ell, $variant_mu, $variant_modt_bits>(
                            rng, &self.0, ct_out,
                        );

                    SharedSecret(shared_secret)
                }
            }

            impl $privkey_name {
                /// Decapsulates an encapsulated key and returns the resulting shared secret. If
                /// the encapsulated key is invalid, then the shared secret will be psuedorandom
                /// garbage.
                pub fn decapsulate(&self, encapsulated_key: &$ciphertext_name) -> SharedSecret {
                    SharedSecret(crate::kem::decap::<
                        $variant_ell,
                        $variant_mu,
                        $variant_modt_bits,
                    >(&self.0, encapsulated_key))
                }
            }

            /// Basic test that keygen, encap, decap, ser, and deser work
            #[test]
            fn test_api() {
                let mut rng = rand::rng();
                let sk = $privkey_name::generate(&mut rng);
                let pk = sk.public_key();

                // Serialize and deserialize the keys
                let mut sk_bytes = [0u8; $privkey_name::SERIALIZED_LEN];
                sk.to_bytes(&mut sk_bytes);
                let sk = $privkey_name::from_bytes(&sk_bytes);

                let mut pk_bytes = [0u8; $pubkey_name::SERIALIZED_LEN];
                pk.to_bytes(&mut pk_bytes);
                let pk = $pubkey_name::from_bytes(&pk_bytes);

                let (ct, ss1) = pk.encapsulate(&mut rng);
                let ct_bytes = ct.as_ref();

                let ct_arr = ct_bytes.try_into().unwrap();
                let ss2 = sk.decapsulate(&ct_arr);

                assert_eq!(ss1.as_bytes(), ss2.as_bytes());
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
    LIGHTSABER_CIPHERTEXT_LEN,
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
    SABER_CIPHERTEXT_LEN,
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
    FIRESABER_CIPHERTEXT_LEN,
    FIRESABER_L,
    FIRESABER_MU,
    FIRESABER_MODULUS_T_BITS
);
