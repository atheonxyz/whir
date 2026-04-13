//! Benchmark comparing power-of-2 (pre-NTT3) vs smooth-{2,3} (NTT3) paths.
//!
//! For each real witness size N = 3 * 2^a:
//!   - pow2 path:   pad N → next_power_of_two(N), commit+prove+verify
//!   - smooth path: commit+prove+verify directly at N (uses ternary folding)
//!
//! Outputs CSV: n,path,prove_ms,prove_peak_mb,verify_ms,verify_peak_mb

use ark_ff::AdditiveGroup;
use ark_std::rand::thread_rng;
use std::alloc::{GlobalAlloc, Layout, System};
use std::borrow::Cow;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;
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

// ─── Tracking allocator ──────────────────────────────────────────────
static CURRENT: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);
struct TA;
unsafe impl GlobalAlloc for TA {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let p = unsafe { System.alloc(l) };
        if !p.is_null() {
            let o = CURRENT.fetch_add(l.size(), Ordering::Relaxed);
            PEAK.fetch_max(o + l.size(), Ordering::Relaxed);
        }
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        unsafe { System.dealloc(p, l) };
        CURRENT.fetch_sub(l.size(), Ordering::Relaxed);
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, ns: usize) -> *mut u8 {
        let np = unsafe { System.realloc(p, l, ns) };
        if !np.is_null() {
            if ns > l.size() {
                let d = ns - l.size();
                let o = CURRENT.fetch_add(d, Ordering::Relaxed);
                PEAK.fetch_max(o + d, Ordering::Relaxed);
            } else {
                CURRENT.fetch_sub(l.size() - ns, Ordering::Relaxed);
            }
        }
        np
    }
}
#[global_allocator]
static ALLOC: TA = TA;
fn reset_peak() {
    PEAK.store(CURRENT.load(Ordering::Relaxed), Ordering::Relaxed);
}
fn get_peak_mb() -> f64 {
    PEAK.load(Ordering::Relaxed) as f64 / (1024.0 * 1024.0)
}

type BF = Field256;
type Cfg = whir::protocols::whir::Config<Identity<BF>>;
const FF: usize = 3;

fn wp() -> ProtocolParameters {
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

/// Run a full prove+verify cycle at `commit_size`, with `real_data` actual witness elements.
/// Returns (prove_ms, prove_peak_mb, verify_ms, verify_peak_mb).
fn run(real_data: usize, commit_size: usize) -> (f64, f64, f64, f64) {
    let p = Cfg::new(commit_size, &wp());
    let mut witness: Vec<BF> = random_vector(thread_rng(), real_data);
    witness.resize(commit_size, BF::ZERO);
    let vectors = vec![witness];
    let vr = vectors.iter().map(|x| x.as_slice()).collect::<Vec<_>>();
    let s = format!("b-{commit_size}");
    let ds = DomainSeparator::protocol(&p).session(&s).instance(&Empty);
    let mut ps = ProverState::new_std(&ds);
    let w = p.commit(&mut ps, &vr);
    let mut lfs: Vec<Box<dyn Evaluate<Identity<BF>>>> = Vec::new();
    lfs.push(Box::new(Covector {
        vector: (0..commit_size as u64).map(BF::from).collect(),
    }));
    let vals = lfs
        .iter()
        .flat_map(|l| vr.iter().map(|v| l.evaluate(p.embedding(), v)))
        .collect::<Vec<_>>();
    let pf: Vec<Box<dyn LinearForm<BF>>> = vec![Box::new(Covector {
        vector: (0..commit_size as u64).map(BF::from).collect(),
    })];

    // ── Prove ──
    reset_peak();
    let t0 = Instant::now();
    let _ = p.prove(
        &mut ps,
        vectors
            .iter()
            .map(|x| Cow::Borrowed(x.as_slice()))
            .collect(),
        vec![Cow::Owned(w)],
        pf,
        Cow::Borrowed(vals.as_slice()),
    );
    let prove_ms = t0.elapsed().as_secs_f64() * 1000.0;
    let prove_peak = get_peak_mb();

    // ── Verify ──
    let proof = ps.proof();
    reset_peak();
    let t1 = Instant::now();
    let mut vs = VerifierState::new_std(&ds, &proof);
    let c = p.receive_commitment(&mut vs).unwrap();
    let wr = lfs
        .iter()
        .map(|w| w.as_ref() as &dyn LinearForm<BF>)
        .collect::<Vec<_>>();
    p.verify(&mut vs, &[&c], &vals).unwrap().verify(wr).unwrap();
    let verify_ms = t1.elapsed().as_secs_f64() * 1000.0;
    let verify_peak = get_peak_mb();

    (prove_ms, prove_peak, verify_ms, verify_peak)
}

fn bench(real_data: usize, commit_size: usize) -> (f64, f64, f64, f64) {
    // Warmup
    let _ = run(real_data, commit_size);
    // Best of 3
    let mut best = run(real_data, commit_size);
    for _ in 0..2 {
        let r = run(real_data, commit_size);
        if r.0 < best.0 {
            best = r;
        }
    }
    best
}

fn next_smooth(mut n: usize) -> usize {
    loop {
        let mut t = n;
        while t % 2 == 0 {
            t /= 2;
        }
        while t % 3 == 0 {
            t /= 3;
        }
        if t == 1 {
            return n;
        }
        n += 1;
    }
}

fn main() {
    let d = 1usize << FF; // interleaving depth must divide size
                          // Target real-data sizes: 3 * 2^a for various a
    let mut targets: Vec<usize> = Vec::new();
    for a in 0..30 {
        let s = 3 * (1 << a);
        if s >= (1 << 12) && s <= (1 << 18) && s % d == 0 {
            targets.push(s);
        }
    }
    targets.sort();
    targets.dedup();

    println!("n,path,prove_ms,prove_peak_mb,verify_ms,verify_peak_mb");
    for &n in &targets {
        // ── Smooth path (NTT3 enabled) ──
        let mut cs = next_smooth(n);
        while cs % d != 0 {
            cs = next_smooth(cs + 1);
        }
        let r = bench(n, cs);
        println!("{n},smooth,{:.3},{:.1},{:.3},{:.1}", r.0, r.1, r.2, r.3);
        eprintln!(
            "  n={n:>7} smooth={cs:>7}  prove={:.1}ms mem={:.0}MB verify={:.2}ms vmem={:.0}MB",
            r.0, r.1, r.2, r.3
        );

        // ── Pow2 path (pre-NTT3 baseline) ──
        let pow2 = n.next_power_of_two();
        let r2 = bench(n, pow2);
        println!("{n},pow2,{:.3},{:.1},{:.3},{:.1}", r2.0, r2.1, r2.2, r2.3);
        eprintln!(
            "  n={n:>7} pow2={pow2:>7}    prove={:.1}ms mem={:.0}MB verify={:.2}ms vmem={:.0}MB",
            r2.0, r2.1, r2.2, r2.3
        );
    }
}
