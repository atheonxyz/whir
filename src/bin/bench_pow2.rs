//! Benchmark for the pow2-only path.
//! For every target data size N, pads to next_power_of_two(N) with zeros.
//! Tracks prove time, prove peak mem, verify time, verify peak mem separately.
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

/// Run prove+verify at `commit_size` (must be pow2) with `real_n` actual values.
/// Returns (prove_ms, prove_peak_mb, verify_ms, verify_peak_mb).
fn run(real_n: usize, commit_size: usize) -> (f64, f64, f64, f64) {
    let p = Cfg::new(commit_size, &wp());
    // Real data padded with zeros
    let mut vec: Vec<BF> = random_vector(thread_rng(), real_n);
    vec.resize(commit_size, BF::ZERO);
    let v = vec![vec];
    let vr = v.iter().map(|x| x.as_slice()).collect::<Vec<_>>();
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
        v.iter().map(|x| Cow::Borrowed(x.as_slice())).collect(),
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

fn bench(real_n: usize, commit_size: usize) -> (f64, f64, f64, f64) {
    let _ = run(real_n, commit_size);
    let mut b = run(real_n, commit_size);
    for _ in 0..2 {
        let r = run(real_n, commit_size);
        if r.0 < b.0 {
            b = r;
        }
    }
    b
}

fn main() {
    let d = 1usize << FF;
    let lo = 1usize << 12;
    let hi = 1usize << 24;

    // Same target sizes as bench_smooth: pow2, 3*2^a, 9*2^a
    let mut targets: Vec<usize> = Vec::new();
    let mut p = lo;
    while p <= hi {
        targets.push(p);
        p *= 2;
    }
    for a in 0..30 {
        let s = 3 * (1 << a);
        if s >= lo && s <= hi && s % d == 0 {
            targets.push(s);
        }
    }
    for a in 0..30 {
        let s = 9 * (1 << a);
        if s >= lo && s <= hi && s % d == 0 {
            targets.push(s);
        }
    }
    targets.sort();
    targets.dedup();

    println!("n,commit_size,prove_ms,prove_peak_mb,verify_ms,verify_peak_mb");
    for &n in &targets {
        let pow2 = n.next_power_of_two();
        let r = bench(n, pow2);
        println!("{n},{pow2},{:.3},{:.1},{:.3},{:.1}", r.0, r.1, r.2, r.3);
        eprintln!(
            "  n={n:>7} -> pow2={pow2:>7}  prove={:.1}ms  pmem={:.0}MB  verify={:.2}ms  vmem={:.0}MB",
            r.0, r.1, r.2, r.3
        );
    }
}
