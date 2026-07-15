//! This module contains code for running known-answer tests (KATs).
//!
//! These tests read test vectors from `test_vectors-kopis<LEVEL>.jsonl` (where `<LEVEL>` is 512,
//! 768, or 1024) and verify that our implementation reproduces the recorded values.
//!
//! Each line of a `.jsonl` file is a JSON object of the following form (all byte strings are
//! encoded as lowercase hex):
//!
//! * `description`:      A text string describing what this test vector is testing.
//! * `sk`:               The KEM secret key.
//! * `pk`:               The serialized KEM public key corresponding to `sk`.
//! * `encap_randomness`: The randomness used to encapsulate to `pk`.
//! * `encapper_ct`:      The serialized KEM ciphertext produced from the encapsulation.
//! * `decapper_ct`:      A KEM ciphertext. This may or may not be equal to `encapper_ct`.
//! * `encapper_ss`:      The KEM shared secret corresponding to the encapsulation that produced
//!                       `encapper_ct`.
//! * `decapper_ss`:      The KEM shared secret from decapsulating `decapper_ct` with `sk`.
//! * `malformed`:        Whether a value above is malformed, making this test vector invalid.
//!
//! A vector with `malformed == true` is skipped, since it does not represent a valid computation.

use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use kopis_kem::{kopis1024, kopis512, kopis768};
use serde::Deserialize;

/// A single known-answer test vector. All byte strings are stored as lowercase hex.
#[derive(Debug, Clone, Deserialize)]
struct KatVector {
    /// A text description of what this test vector exercises.
    description: String,
    /// The KEM secret key (32-byte seed).
    #[serde(with = "hex::serde")]
    sk: Vec<u8>,
    /// The serialized KEM public key corresponding to `sk`.
    #[serde(with = "hex::serde")]
    pk: Vec<u8>,
    /// The randomness used to encapsulate to `pk`.
    #[serde(with = "hex::serde")]
    encap_randomness: Vec<u8>,
    /// The serialized ciphertext produced by the encapsulation.
    #[serde(with = "hex::serde")]
    encapper_ct: Vec<u8>,
    /// The ciphertext handed to the decapsulator. May differ from `encapper_ct`.
    #[serde(with = "hex::serde")]
    decapper_ct: Vec<u8>,
    /// The shared secret from the encapsulation that produced `encapper_ct`.
    #[serde(with = "hex::serde")]
    encapper_ss: Vec<u8>,
    /// The shared secret from decapsulating `decapper_ct` with `sk`.
    #[serde(with = "hex::serde")]
    decapper_ss: Vec<u8>,
    /// Whether a value above is malformed, making this test vector invalid.
    malformed: bool,
}

/// Reads vectors from `path`, one JSON object per line.
fn read_vectors(path: &Path) -> Vec<KatVector> {
    let file = File::open(path)
        .unwrap_or_else(|e| panic!("could not open test-vector file {}: {e}", path.display()));
    BufReader::new(file)
        .lines()
        .enumerate()
        .filter_map(|(i, line)| {
            let line = line.expect("could not read test-vector line");
            // Skip blank lines so trailing newlines don't cause parse errors.
            if line.trim().is_empty() {
                return None;
            }
            let v = serde_json::from_str(&line)
                .unwrap_or_else(|e| panic!("could not parse vector on line {}: {e}", i + 1));
            Some(v)
        })
        .collect()
}

/// Reads and verifies every (non-malformed) vector for a single Kopis level. This is a macro
/// because each level uses distinct key/ciphertext types.
macro_rules! kat_test {
    (
        $test_name:ident,
        $level:expr,
        $sk_ty:ty,
        $pk_ty:ty,
        $ct_len:expr
    ) => {
        #[test]
        fn $test_name() {
            let path_str = format!("tests/test_vectors-kopis{}.jsonl", $level);
            let path = Path::new(&path_str);
            let vectors = read_vectors(path);

            for vector in &vectors {
                if vector.malformed {
                    continue;
                }

                let ctx = || format!("Kopis-{} vector failed: {}", $level, vector.description);

                // Expand the secret key from its 32-byte seed and check the derived public key
                // matches the recorded one.
                let seed: [u8; 32] = vector
                    .sk
                    .as_slice()
                    .try_into()
                    .unwrap_or_else(|_| panic!("{}: sk is not 32 bytes", ctx()));
                let sk = <$sk_ty>::expand_from_seed(&seed);
                let pk = sk.public_key();

                let mut pk_bytes = [0u8; <$pk_ty>::SERIALIZED_LEN];
                pk.serialize(&mut pk_bytes);
                assert_eq!(
                    pk_bytes.as_slice(),
                    vector.pk.as_slice(),
                    "{}: derived public key does not match recorded pk",
                    ctx()
                );

                // Re-run the (deterministic) encapsulation and check the ciphertext and shared
                // secret match the recorded encapper values.
                let encap_randomness: [u8; 32] = vector
                    .encap_randomness
                    .as_slice()
                    .try_into()
                    .unwrap_or_else(|_| panic!("{}: encap_randomness is not 32 bytes", ctx()));
                let (encapper_ct, encapper_ss) = pk.encapsulate_deterministic(&encap_randomness);
                assert_eq!(
                    encapper_ct.as_slice(),
                    vector.encapper_ct.as_slice(),
                    "{}: recomputed encapper_ct does not match recorded value",
                    ctx()
                );
                assert_eq!(
                    encapper_ss.as_bytes().as_slice(),
                    vector.encapper_ss.as_slice(),
                    "{}: recomputed encapper_ss does not match recorded value",
                    ctx()
                );

                // Decapsulate the recorded decapper ciphertext and check the shared secret matches.
                let decapper_ct: &[u8; $ct_len] = vector
                    .decapper_ct
                    .as_slice()
                    .try_into()
                    .unwrap_or_else(|_| panic!("{}: decapper_ct has wrong length", ctx()));
                let decapper_ss = sk.decapsulate(decapper_ct);
                assert_eq!(
                    decapper_ss.as_bytes().as_slice(),
                    vector.decapper_ss.as_slice(),
                    "{}: recomputed decapper_ss does not match recorded value",
                    ctx()
                );
            }
        }
    };
}

kat_test!(
    kat_kopis512,
    512,
    kopis512::Kopis512SecretKey,
    kopis512::Kopis512PublicKey,
    kopis512::KOPIS512_CIPHERTEXT_LEN
);
kat_test!(
    kat_kopis768,
    768,
    kopis768::Kopis768SecretKey,
    kopis768::Kopis768PublicKey,
    kopis768::KOPIS768_CIPHERTEXT_LEN
);
kat_test!(
    kat_kopis1024,
    1024,
    kopis1024::Kopis1024SecretKey,
    kopis1024::Kopis1024PublicKey,
    kopis1024::KOPIS1024_CIPHERTEXT_LEN
);
