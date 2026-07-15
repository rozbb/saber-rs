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

/// Defines convenience types and impls for a given Kopis variant
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
                /// Generate a fresh secret key
                pub fn generate(rng: &mut impl CryptoRng) -> Self {
                    Self(KemSecretKey::generate::<$variant_mu>(rng))
                }

                /// Returns the seed that produced this secret key
                pub fn seed(&self) -> [u8; 32] {
                    self.0.seed()
                }

                /// Deserializes a secret key from a 32-byte seed
                pub fn expand_from_seed(bytes: &[u8; 32]) -> Self {
                    Self(KemSecretKey::expand_from_seed::<$variant_mu>(bytes))
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
                pub fn serialize(&self, out_buf: &mut [u8; Self::SERIALIZED_LEN]) {
                    self.0.serialize(out_buf);
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
                    let ss = crate::kem::encap::<$variant_ell, $variant_mu, $variant_modt_bits>(
                        rng, &self.0, &mut ct,
                    );

                    (ct, SharedSecret(ss))
                }
            }

            impl $privkey_name {
                /// Decapsulates an encapsulated key and returns the resulting shared secret. If
                /// the encapsulated key is invalid, then the shared secret will be pseudorandom
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
                let sk_seed = sk.seed();
                let sk = $privkey_name::expand_from_seed(&sk_seed);

                let mut pk_bytes = [0u8; $pubkey_name::SERIALIZED_LEN];
                pk.serialize(&mut pk_bytes);
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
    kopis512,
    "Kopis-512 is designed to have security close to that of AES-128",
    Kopis512PublicKey,
    Kopis512SecretKey,
    Kopis512Ciphertext,
    KOPIS512_CIPHERTEXT_LEN,
    KOPIS512_L,
    KOPIS512_MU,
    KOPIS512_T
);

variant_impl!(
    kopis768,
    "Kopis-768 is designed to have security close to that of AES-192",
    Kopis768PublicKey,
    Kopis768SecretKey,
    Kopis768Ciphertext,
    KOPIS768_CIPHERTEXT_LEN,
    KOPIS768_L,
    KOPIS768_MU,
    KOPIS768_T
);

variant_impl!(
    kopis1024,
    "Kopis-1024 is designed to have security close to that of AES-256",
    Kopis1024PublicKey,
    Kopis1024SecretKey,
    Kopis1024Ciphertext,
    KOPIS1024_CIPHERTEXT_LEN,
    KOPIS1024_L,
    KOPIS1024_MU,
    KOPIS1024_T
);
