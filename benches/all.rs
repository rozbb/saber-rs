use saber_kem::{
    firesaber::{FiresaberPublicKey, FiresaberSecretKey},
    lightsaber::{LightsaberPublicKey, LightsaberSecretKey},
    saber::{SaberPublicKey, SaberSecretKey},
};

use criterion::{criterion_group, criterion_main, Criterion};

macro_rules! bench_variant {
    ($bench_name:ident, $privkey_name:ident, $pubkey_name:ident) => {
        fn $bench_name(c: &mut Criterion) {
            let mut rng = rand::rng();

            let gen_bench_name = format!("{}-gen-keypair", stringify!($bench_name));
            c.bench_function(&gen_bench_name, |b| {
                b.iter(|| $privkey_name::generate(&mut rng))
            });
            let sk = $privkey_name::generate(&mut rng);
            let pk = sk.public_key();
            let pk_bytes = {
                let mut buf = [0u8; $pubkey_name::SERIALIZED_LEN];
                pk.to_bytes(&mut buf);
                buf
            };

            let encap_bench_name = format!("{}-deser-pubkey", stringify!($bench_name));
            c.bench_function(&encap_bench_name, |b| {
                b.iter(|| $pubkey_name::from_bytes(&pk_bytes))
            });
            let (ct, _) = pk.encapsulate(&mut rng);

            let encap_bench_name = format!("{}-encap", stringify!($bench_name));
            c.bench_function(&encap_bench_name, |b| b.iter(|| pk.encapsulate(&mut rng)));
            let (ct, _) = pk.encapsulate(&mut rng);

            let decap_bench_name = format!("{}-decap", stringify!($bench_name));
            c.bench_function(&decap_bench_name, |b| b.iter(|| sk.decapsulate(&ct)));
        }
    };
}

bench_variant!(lightsaber, LightsaberSecretKey, LightsaberPublicKey);
bench_variant!(saber, SaberSecretKey, SaberPublicKey);
bench_variant!(firesaber, FiresaberSecretKey, FiresaberPublicKey);

criterion_group!(benches, lightsaber, saber, firesaber);
criterion_main!(benches);
