/// Simple profiling harness: runs Saber KEM keygen/encap/decap in a tight loop.
/// Use with `cargo flamegraph --example profile` or `sample` to find hot spots.
fn main() {
    let mut rng = rand::rng();
    let iterations = 10_000;

    for _ in 0..iterations {
        let sk = saber_kem::saber::SaberSecretKey::generate(&mut rng);
        let pk = sk.public_key();
        let (ct, _ss1) = pk.encapsulate(&mut rng);
        let _ss2 = sk.decapsulate(&ct);
    }
}
