//! Sweep 2^a * 3^b (a<=23, b<=2): smooth+ternary vs pow2-padded.
//!
//! Variant A (pow2): commit next_power_of_two(N) with zero-pad.
//! Variant B (smooth_tern): commit N directly, mixed radix-2/3 folding.
//!
//! Writes CSV incrementally so a crash/OOM on the tail does not lose earlier rows.

use ark_ff::AdditiveGroup;
use ark_std::rand::thread_rng;
use std::alloc::{GlobalAlloc, Layout, System};
use std::borrow::Cow;
use std::fs::OpenOptions;
use std::io::Write;
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

// ---- Allocation tracker ----
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
fn peak_mb() -> f64 {
    PEAK.load(Ordering::Relaxed) as f64 / (1024.0 * 1024.0)
}

type BF = Field256;
type Cfg = whir::protocols::whir::Config<Identity<BF>>;
const FF: usize = 3;

fn wp() -> ProtocolParameters {
    ProtocolParameters {
        security_level: 128,
        pow_bits: 10,
        initial_folding_factor: FF,
        folding_factor: FF,
        unique_decoding: false,
        starting_log_inv_rate: 2,
        batch_size: 1,
        hash_id: hash::SHA2,
    }
}

struct Res {
    prove_ms: f64,
    prove_mb: f64,
    verify_ms: f64,
    verify_mb: f64,
}

fn run_once(n: usize, commit_size: usize) -> Res {
    let cfg = Cfg::new(commit_size, &wp());

    let mut witness: Vec<BF> = random_vector(thread_rng(), n);
    witness.resize(commit_size, BF::ZERO);
    let vectors = vec![witness];
    let vec_refs = vectors.iter().map(|v| v.as_slice()).collect::<Vec<_>>();
    let session = format!("b-{commit_size}");
    let ds = DomainSeparator::protocol(&cfg).session(&session).instance(&Empty);

    let mut ps = ProverState::new_std(&ds);
    let w = cfg.commit(&mut ps, &vec_refs);
    let lfs: Vec<Box<dyn Evaluate<Identity<BF>>>> = vec![Box::new(Covector {
        vector: (0..commit_size as u64).map(BF::from).collect(),
    })];
    let vals = lfs
        .iter()
        .flat_map(|l| vec_refs.iter().map(|v| l.evaluate(cfg.embedding(), v)))
        .collect::<Vec<_>>();
    let pf: Vec<Box<dyn LinearForm<BF>>> = vec![Box::new(Covector {
        vector: (0..commit_size as u64).map(BF::from).collect(),
    })];

    reset_peak();
    let t0 = Instant::now();
    let _ = cfg.prove(
        &mut ps,
        vectors.iter().map(|v| Cow::Borrowed(v.as_slice())).collect(),
        vec![Cow::Owned(w)],
        pf,
        Cow::Borrowed(vals.as_slice()),
    );
    let prove_ms = t0.elapsed().as_secs_f64() * 1000.0;
    let prove_mb = peak_mb();

    let proof = ps.proof();
    reset_peak();
    let t1 = Instant::now();
    let mut vs = VerifierState::new_std(&ds, &proof);
    let c = cfg.receive_commitment(&mut vs).unwrap();
    let wr = lfs.iter().map(|w| w.as_ref() as &dyn LinearForm<BF>).collect::<Vec<_>>();
    cfg.verify(&mut vs, &[&c], &vals).unwrap().verify(wr).unwrap();
    let verify_ms = t1.elapsed().as_secs_f64() * 1000.0;
    let verify_mb = peak_mb();

    Res { prove_ms, prove_mb, verify_ms, verify_mb }
}

fn best_of(n: usize, commit_size: usize, iters: u32) -> Res {
    let _ = run_once(n, commit_size); // warmup
    let mut b = run_once(n, commit_size);
    for _ in 0..iters.saturating_sub(1) {
        let r = run_once(n, commit_size);
        if r.prove_ms + r.verify_ms < b.prove_ms + b.verify_ms {
            b = r;
        }
    }
    b
}

fn next_pow2(n: usize) -> usize {
    if n.is_power_of_two() { n } else { n.next_power_of_two() }
}

fn main() {
    // Env-controlled caps so tail can be skipped if OOM / too slow.
    let a_max: u32 = std::env::var("A_MAX").ok().and_then(|s| s.parse().ok()).unwrap_or(23);
    let b_max: u32 = std::env::var("B_MAX").ok().and_then(|s| s.parse().ok()).unwrap_or(2);
    let csv_path = std::env::var("CSV").unwrap_or_else(|_| "test-arth/bench_smooth_vs_pow2_sweep.csv".into());
    let iters_small: u32 = std::env::var("ITERS_SMALL").ok().and_then(|s| s.parse().ok()).unwrap_or(3);
    let iters_large: u32 = std::env::var("ITERS_LARGE").ok().and_then(|s| s.parse().ok()).unwrap_or(1);
    let large_cutoff: usize = 1 << 20;

    // Build size list: 2^a · 3^b with a≥3 (need 8 | size) and in ranges requested.
    let mut sizes: Vec<(usize, u32, u32)> = Vec::new();
    for b in 0..=b_max {
        let three_pow = 3u64.pow(b) as usize;
        for a in 3..=a_max {
            let n = three_pow.checked_shl(a).unwrap_or(0);
            if n == 0 { continue; }
            if n % (1 << FF) != 0 { continue; }
            sizes.push((n, a, b));
        }
    }
    sizes.sort_by_key(|x| x.0);
    sizes.dedup_by_key(|x| x.0);

    let mut f = OpenOptions::new().create(true).truncate(true).write(true).open(&csv_path).unwrap();
    writeln!(f, "n,a,b,variant,commit_size,prove_ms,prove_mb,verify_ms,verify_mb").unwrap();
    println!("n,a,b,variant,commit_size,prove_ms,prove_mb,verify_ms,verify_mb");
    eprintln!(
        "{:>12} {:>3} {:>2} | {:>12} | {:>10} | {:>9} {:>7} {:>9} {:>7}",
        "n","a","b","variant","commit_sz","prove_ms","pmb","verify_ms","vmb"
    );

    for &(n, a, b) in &sizes {
        let iters = if n >= large_cutoff { iters_large } else { iters_small };
        let pow2_n = next_pow2(n);

        // A) pow2
        let r = best_of(n, pow2_n, iters);
        let line = format!("{n},{a},{b},pow2,{pow2_n},{:.3},{:.2},{:.3},{:.2}",
            r.prove_ms, r.prove_mb, r.verify_ms, r.verify_mb);
        println!("{line}");
        writeln!(f, "{line}").unwrap();
        f.flush().unwrap();
        eprintln!(
            "{:>12} {:>3} {:>2} | {:>12} | {:>10} | {:>9.1} {:>7.1} {:>9.3} {:>7.1}",
            n,a,b,"pow2",pow2_n,r.prove_ms,r.prove_mb,r.verify_ms,r.verify_mb
        );

        // B) smooth + ternary
        let r = best_of(n, n, iters);
        let line = format!("{n},{a},{b},smooth_tern,{n},{:.3},{:.2},{:.3},{:.2}",
            r.prove_ms, r.prove_mb, r.verify_ms, r.verify_mb);
        println!("{line}");
        writeln!(f, "{line}").unwrap();
        f.flush().unwrap();
        eprintln!(
            "{:>12} {:>3} {:>2} | {:>12} | {:>10} | {:>9.1} {:>7.1} {:>9.3} {:>7.1}",
            n,a,b,"smooth_tern",n,r.prove_ms,r.prove_mb,r.verify_ms,r.verify_mb
        );
    }
}
