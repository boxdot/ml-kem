//! ML-KEM-768 round-trip benchmark.
//!
//! Source <https://github.com/FiloSottile/mlkem768/blob/main/mlkem768_test.go> BenchmarkRoundTrip.

use criterion::{Criterion, black_box, criterion_group, criterion_main};
use ml_kem_example::kem::mlkem::EncapsulationKey;
use sha3::{
    Shake128,
    digest::{ExtendableOutput, Update, XofReader},
};

/// Deterministic 32-byte stream so each iteration uses fresh-looking randomness without any
/// syscalls or crypto-rng dependency.
struct Rng(sha3::Shake128Reader);

impl Rng {
    fn new(label: &[u8]) -> Self {
        let mut x = Shake128::default();
        x.update(label);
        Self(x.finalize_xof())
    }

    fn bytes32(&mut self) -> [u8; 32] {
        let mut out = [0u8; 32];
        self.0.read(&mut out);
        out
    }
}

fn bench_roundtrip(c: &mut Criterion) {
    // Fixed long-lived keypair (Alice's)
    let mut setup = Rng::new(b"mlkem768-bench-setup");
    let d = setup.bytes32();
    let z = setup.bytes32();
    let (ek, dk) = EncapsulationKey::generate(&d, &z);

    // Fixed ciphertext against that keypair (used by the Alice/Decapsulate arm).
    let m = setup.bytes32();
    let (_k, ct) = ek.encaps(&m);

    // Alice: generate a new keypair + decapsulate a known ciphertext.
    c.bench_function("roundtrip/alice", |b| {
        let mut rng = Rng::new(b"mlkem768-bench-alice");
        b.iter(|| {
            let d = rng.bytes32();
            let z = rng.bytes32();
            let (ek_s, _dk_s) = EncapsulationKey::generate(&d, &z);
            black_box(&ek_s);

            let ks = dk.decaps(black_box(&ct));
            black_box(ks);
        });
    });

    // Bob: encapsulate to Alice's public key.
    c.bench_function("roundtrip/bob", |b| {
        let mut rng = Rng::new(b"mlkem768-bench-bob");
        b.iter(|| {
            let m = rng.bytes32();
            let (ks, cs) = ek.encaps(black_box(&m));
            black_box((ks, cs));
        });
    });
}

criterion_group!(benches, bench_roundtrip);
criterion_main!(benches);
