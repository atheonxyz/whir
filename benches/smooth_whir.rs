//! Benchmark: smooth-{2,3} sizes vs power-of-2 padding.
//!
//! For a fixed amount of real witness data (N elements), compare:
//!   - smooth path:  commit N elements directly (N is smooth)
//!   - pow2 path:    pad N to next_power_of_two(N), commit padded
//!
//! Uses BN254 (Field256) which has 3^2 | (p-1).

use std::borrow::Cow;

use ark_ff::AdditiveGroup;
use ark_std::rand::thread_rng;
use divan::{black_box, AllocProfiler, Bencher};

use whir::{
    algebra::{
        embedding::Identity,
        fields::Field256,
        linear_form::{Covector, Evaluate, LinearForm},
        random_vector,
    },
    hash,
    parameters::ProtocolParameters,
    transcript::{codecs::Empty, DomainSeparator, ProverState, VerifierState},
};

#[global_allocator]
static ALLOC: AllocProfiler = AllocProfiler::system();

type BF = Field256;
type Cfg = whir::protocols::whir::Config<Identity<BF>>;

const FF: usize = 3;

fn whir_params() -> ProtocolParameters {
    ProtocolParameters {
        security_level: 32,
        pow_bits: 0,
        initial_folding_factor: FF,
        folding_factor: FF,
        unique_decoding: false,
        starting_log_inv_rate: 1,
        batch_size: 1,
        hash_id: hash::SHA2,
    }
}

/// Run full WHIR prove+verify.
/// `real_data` = number of actual witness elements.
/// `commit_size` = the size passed to Config (smooth or pow2), >= real_data.
/// Extra positions are zero-padded.
fn prove_and_verify(real_data: usize, commit_size: usize) {
    let params = Cfg::new(commit_size, &whir_params());

    // Real witness data + zero padding
    let mut witness: Vec<BF> = random_vector(thread_rng(), real_data);
    witness.resize(commit_size, BF::ZERO);
    let vectors = vec![witness];
    let vec_refs = vectors.iter().map(|v| v.as_slice()).collect::<Vec<_>>();

    let session = format!("bench-{commit_size}");
    let ds = DomainSeparator::protocol(&params)
        .session(&session)
        .instance(&Empty);

    let mut prover_state = ProverState::new_std(&ds);
    let w = params.commit(&mut prover_state, &vec_refs);

    let mut lfs: Vec<Box<dyn Evaluate<Identity<BF>>>> = Vec::new();
    lfs.push(Box::new(Covector {
        vector: (0..commit_size as u64).map(BF::from).collect(),
    }));
    let values = lfs
        .iter()
        .flat_map(|lf| vec_refs.iter().map(|v| lf.evaluate(params.embedding(), v)))
        .collect::<Vec<_>>();

    let prove_forms: Vec<Box<dyn LinearForm<BF>>> = vec![Box::new(Covector {
        vector: (0..commit_size as u64).map(BF::from).collect(),
    })];

    let _ = params.prove(
        &mut prover_state,
        vectors
            .iter()
            .map(|v| Cow::Borrowed(v.as_slice()))
            .collect(),
        vec![Cow::Owned(w)],
        prove_forms,
        Cow::Borrowed(values.as_slice()),
    );

    let proof = prover_state.proof();
    let mut vs = VerifierState::new_std(&ds, &proof);
    let c = params.receive_commitment(&mut vs).unwrap();
    let w_refs = lfs
        .iter()
        .map(|w| w.as_ref() as &dyn LinearForm<BF>)
        .collect::<Vec<_>>();
    params
        .verify(&mut vs, &[&c], &values)
        .unwrap()
        .verify(w_refs)
        .unwrap();
}

// ─── Each pair: same real data, smooth commit vs pow2 commit ────────
//
// For N real elements where N = 3 * 2^a:
//   smooth commits N elements (no waste)
//   pow2 commits next_power_of_two(N) = 2^(a+2) elements (33% waste)
//
// The smooth path does less NTT, less Merkle, less sumcheck work.

// --- N = 768 = 3*256 (pow2 = 1024, 25% waste) ---

#[divan::bench]
fn n768_smooth_768(bencher: Bencher) {
    bencher.bench(|| black_box(prove_and_verify(768, 768)));
}
#[divan::bench]
fn n768_pow2_1024(bencher: Bencher) {
    bencher.bench(|| black_box(prove_and_verify(768, 1024)));
}

// --- N = 1536 = 3*512 (pow2 = 2048, 25% waste) ---

#[divan::bench]
fn n1536_smooth_1536(bencher: Bencher) {
    bencher.bench(|| black_box(prove_and_verify(1536, 1536)));
}
#[divan::bench]
fn n1536_pow2_2048(bencher: Bencher) {
    bencher.bench(|| black_box(prove_and_verify(1536, 2048)));
}

// --- N = 3072 = 3*1024 (pow2 = 4096, 25% waste) ---

#[divan::bench]
fn n3072_smooth_3072(bencher: Bencher) {
    bencher.bench(|| black_box(prove_and_verify(3072, 3072)));
}
#[divan::bench]
fn n3072_pow2_4096(bencher: Bencher) {
    bencher.bench(|| black_box(prove_and_verify(3072, 4096)));
}

// --- N = 6144 = 3*2048 (pow2 = 8192, 25% waste) ---

#[divan::bench]
fn n6144_smooth_6144(bencher: Bencher) {
    bencher.bench(|| black_box(prove_and_verify(6144, 6144)));
}
#[divan::bench]
fn n6144_pow2_8192(bencher: Bencher) {
    bencher.bench(|| black_box(prove_and_verify(6144, 8192)));
}

// --- N = 12288 = 3*4096 (pow2 = 16384, 25% waste) ---

#[divan::bench]
fn n12288_smooth_12288(bencher: Bencher) {
    bencher.bench(|| black_box(prove_and_verify(12288, 12288)));
}
#[divan::bench]
fn n12288_pow2_16384(bencher: Bencher) {
    bencher.bench(|| black_box(prove_and_verify(12288, 16384)));
}

fn main() {
    divan::main();
}
