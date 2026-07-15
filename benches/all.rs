use kopis_kem::{
    kopis1024::Kopis1024SecretKey, kopis512::Kopis512SecretKey, kopis768::Kopis768SecretKey,
};

use criterion::{criterion_group, criterion_main, Criterion};

macro_rules! bench_variant {
    ($bench_name:ident, $privkey_name:ident) => {
        fn $bench_name(c: &mut Criterion) {
            let mut rng = rand::rng();

            let gen_bench_name = format!("{}-gen-keypair", stringify!($bench_name));
            c.bench_function(&gen_bench_name, |b| {
                b.iter(|| $privkey_name::generate(&mut rng))
            });
            let sk = $privkey_name::generate(&mut rng);
            let pk = sk.public_key();

            let encap_bench_name = format!("{}-encap", stringify!($bench_name));
            c.bench_function(&encap_bench_name, |b| b.iter(|| pk.encapsulate(&mut rng)));
            let (ct, _) = pk.encapsulate(&mut rng);

            let decap_bench_name = format!("{}-decap", stringify!($bench_name));
            c.bench_function(&decap_bench_name, |b| b.iter(|| sk.decapsulate(&ct)));
        }
    };
}

bench_variant!(kopis512, Kopis512SecretKey);
bench_variant!(kopis768, Kopis768SecretKey);
bench_variant!(kopis1024, Kopis1024SecretKey);

criterion_group!(benches, kopis512, kopis768, kopis1024);
criterion_main!(benches);
