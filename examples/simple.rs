use saber_kem::{
    kem_traits::{Decapsulate, Encapsulate},
    lightsaber::{LightsaberCiphertext, LightsaberPublicKey, LightsaberSecretKey},
};

fn main() {
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
}
