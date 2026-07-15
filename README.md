kopis_kem
=========

This crate is a pure-Rust, no-std implementation of the Kopis key encapsulation mechanism (KEM). Kopis is a lattice-based KEM that is designed to be secure against classical and quantum adversaries. It comes in three variants:

* Kopis-512, which is designed to have security roughly equivalent to AES-128
* Kopis-768, which is designed to have security roughly equivalent to AES-192
* Kopis-1024, which is designed to have security roughly equivalent to AES-256

Warning
-------

This crate has not been audited in any sense of the word. Use at your own risk.

Example code
------------

The following code can be found in [`examples/simple.rs`](examples/simple.rs).

```rust
use kopis_kem::kopis512::{
    Kopis512Ciphertext, Kopis512PublicKey, Kopis512SecretKey, KOPIS512_CIPHERTEXT_LEN,
};

let mut rng = rand::rng();

// Generate a keypair
let sk = Kopis512SecretKey::generate(&mut rng);
let pk = sk.public_key();

// Serialize the secret key, maybe to save on disk
let sk_seed: [u8; 32] = sk.seed();

// Deserialize the secret key
let sk = Kopis512SecretKey::expand_from_seed(&sk_seed);

// Also serialize and deserialize the public key
let mut pk_bytes = [0u8; Kopis512PublicKey::SERIALIZED_LEN];
pk.to_bytes(&mut pk_bytes);
let slice_containing_pk = pk_bytes.as_slice();
// The API only accepts fixed-len slices, so we have to cast it first
assert_eq!(
    slice_containing_pk.len(),
    Kopis512PublicKey::SERIALIZED_LEN
);
let pk_arr = slice_containing_pk.try_into().unwrap();
let pk = Kopis512PublicKey::from_bytes(pk_arr);

// Encapsulate a shared secret, ss1, to pk
let (_ct, _ss1) = pk.encapsulate(&mut rng);
// Alternatively, if you have a buffer and want to avoid an extra allocation, encapsulate in
// place. Kopis512Ciphertext is just a byte array, so no conversion necessary:
let mut ct = [0u8; KOPIS512_CIPHERTEXT_LEN];
let ss1 = pk.encapsulate_in_place(&mut rng, &mut ct);
let slice_containing_ct = ct.as_slice();

// Deserializing is also straightforward
assert_eq!(slice_containing_ct.len(), KOPIS512_CIPHERTEXT_LEN);
let receiver_ct: &Kopis512Ciphertext = slice_containing_ct.try_into().unwrap();

// Use the secret key to decapsulate the ciphertext
let ss2 = sk.decapsulate(receiver_ct);

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
