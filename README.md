saber_kem
=========

This crate is a pure-Rust, no-std implementation of [draft 3](https://www.esat.kuleuven.be/cosic/pqcrypto/saber/files/saberspecround3.pdf) of the Saber key encapsulation mechanism (KEM). Saber is a lattice-based KEM that is designed to be secure against classical and quantum adversaries. It comes three variants:

* LightSaber, which is designed to have security roughly equivalent to AES-128
* Saber, which is designed to have security roughly equivalent to AES-192
* FireSaber, which is designed to have security roughly equivalent to AES-256

Warning
-------

This crate has not been audited in any sense of the word. Use at your own risk.

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

// Deserialize the secret key
let sk = LightsaberSecretKey::from_bytes(&sk_bytes);

// Also serialize and deserialize the public key
let mut pk_bytes = [0u8; LightsaberPublicKey::SERIALIZED_LEN];
pk.to_bytes(&mut pk_bytes);
let pk = LightsaberPublicKey::from_bytes(&pk_bytes);

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

// Ensure the shared secrets are equal
assert_eq!(ss1, ss2);

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
