use criterion::{criterion_group, criterion_main, Criterion};
use rand::thread_rng;
use saber::RingElem;

pub fn mul(c: &mut Criterion) {
    let mut rng = thread_rng();
    c.bench_function("schoolbook-mul", |b| {
        b.iter(|| {
            let r1 = RingElem::rand(&mut rng);
            let r2 = RingElem::rand(&mut rng);
            r1.schoolbook_mul(&r2)
        })
    });
    c.bench_function("karatsuba-mul", |b| {
        b.iter(|| {
            let r1 = RingElem::rand(&mut rng);
            let r2 = RingElem::rand(&mut rng);
            r1.karatsuba_mul(&r2)
        })
    });
}

criterion_group!(benches, mul);
criterion_main!(benches);