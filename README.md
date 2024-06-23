saber_kem
=========

This crate is a pure-Rust, no-std implementation of [draft 3](https://www.esat.kuleuven.be/cosic/pqcrypto/saber/files/saberspecround3.pdf) of the Saber key encapsulation mechanism (KEM). Saber is a lattice-based KEM that is designed to be secure against classical and quantum adversaries. It comes three variants:

* LightSaber, which is designed to have security roughly equivalent to AES-128
* Saber, which is designed to have security roughly equivalent to AES-192
* FireSaber, which is designed to have security roughly equivalent to AES-256

Warning
-------

This crate has not been audited in any sense of the word. Use at your own risk.

Why Saber?
----------

In general, if you are looking to use a post-quantum KEM and have no other requirements, you should use ML-KEM (aka "Kyber", its pre-standardization name), since it is faster and more standardized than Saber. However, Saber has two small benefits over Kyber:

* All Saber public keys and ciphertexts pack perfectly into bytes. So if you need to perform a keyed permutation on a KEM's public keys, as is required in some ideal-cipher-based constructions such as [CAKE](https://eprint.iacr.org/2023/470), you can simply use a wide-block cipher over a serialized Saber public key. In comparison Kyber requires you to define a permutation over arrays of mod-q values (note: Kyber public keys actually can be compressed to pack into bytes, but nobody has proven it secure; Theorem 2 of the [original paper](https://eprint.iacr.org/2017/634) only considers the uncompressed scheme).
* Relatedly, all Saber arithmetic is modulo a power of two, which is extremely simple for CPUs to work with. Arithmetic modulo a prime can yield much faster computations, but it can also cause accidental timing leaks due to compilers being too smart. Such vulnerabilities have affected [Kyber](https://groups.google.com/a/list.nist.gov/g/pqc-forum/c/hqbtIGFKIpU/m/cnE3pbueBgAJ) and [curve25519](https://rustsec.org/advisories/RUSTSEC-2024-0344.html). I don't claim any of these other projects are insecure, just that this is a specific issue they must contend with going forward, that Saber does not have to.

Compatibility
-------------

This crate is compatible with Saber's [C reference implementation](https://github.com/KULeuven-COSIC/SABER/tree/f7f39e4db2f3e22a21e1dd635e0601caae2b4510). Known-answer tests (KATs) test vectors can be found in [`tests/`](tests/). Test vectors were taken directly from the previously linked repo, and converted to JSON using [`tests/convert_rsp_to_json.py`](tests/convert_rsp_to_json.py).

Example code
------------

The following code can be found in [`examples/simple.rs`](examples/simple.rs).

```rust
use saber_kem::{
    kem_traits::{Decapsulate, Encapsulate},
    lightsaber::{LightsaberCiphertext, LightsaberPublicKey, LightsaberSecretKey},
};

let mut rng = rand::thread_rng();

// Generate a keypair
let sk = LightsaberSecretKey::generate(&mut rng);
let pk = sk.public_key();

// Serialize the secret key, maybe to save on disk
let mut sk_bytes = [0u8; LightsaberSecretKey::SERIALIZED_LEN];
sk.to_bytes(&mut sk_bytes);
let slice_containing_sk = sk_bytes.as_slice();

// Deserialize the secret key
// The API only accepts fixed-len slices, so we have to cast it first
let sk_arr = slice_containing_sk[..LightsaberSecretKey::SERIALIZED_LEN]
    .try_into()
    .unwrap();
let sk = LightsaberSecretKey::from_bytes(sk_arr);

// Also serialize and deserialize the public key
let mut pk_bytes = [0u8; LightsaberPublicKey::SERIALIZED_LEN];
pk.to_bytes(&mut pk_bytes);
let slice_containing_pk = pk_bytes.as_slice();
// The API only accepts fixed-len slices, so we have to cast it first
let pk_arr = slice_containing_pk[..LightsaberPublicKey::SERIALIZED_LEN]
    .try_into()
    .unwrap();
let pk = LightsaberPublicKey::from_bytes(pk_arr);

// Encapsulate a shared secret, ss1, to pk
let (ct, ss1) = pk.encapsulate(&mut rng).unwrap();
// The ciphertext is just bytes, so serializing is straightforward
let ct_bytes = ct.as_ref();

// Deserializing is also straightforward
assert_eq!(ct_bytes.len(), LightsaberCiphertext::LEN);
let mut receiver_ct = LightsaberCiphertext::default();
receiver_ct.as_mut().copy_from_slice(ct_bytes);

// Use the secret key to decapsulate the ciphertext
let ss2 = sk.decapsulate(&receiver_ct).unwrap();

// Check the shared secrets are equal. NOTE is not a constant-time check (ie not secure). We
// only do this for testing purposes.
assert_eq!(ss1.as_bytes(), ss2.as_bytes());

println!("KEM ran successfully");
```

Benchmarks
----------

We have implemented benchmarks for key generation, encapsulation, and decapsulation for all variants. Simply run `cargo bench`.

License
-------

Licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE))
 * MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.
