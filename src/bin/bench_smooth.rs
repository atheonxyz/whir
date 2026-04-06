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

fn run(sz: usize) -> (f64, f64, f64) {
    let p = Cfg::new(sz, &wp());
    let v: Vec<Vec<BF>> = vec![random_vector(thread_rng(), sz)];
    let vr = v.iter().map(|x| x.as_slice()).collect::<Vec<_>>();
    let s = format!("b-{sz}");
    let ds = DomainSeparator::protocol(&p).session(&s).instance(&Empty);
    let mut ps = ProverState::new_std(&ds);
    let w = p.commit(&mut ps, &vr);
    let mut lfs: Vec<Box<dyn Evaluate<Identity<BF>>>> = Vec::new();
    lfs.push(Box::new(Covector {
        vector: (0..sz as u64).map(BF::from).collect(),
    }));
    let vals = lfs
        .iter()
        .flat_map(|l| vr.iter().map(|v| l.evaluate(p.embedding(), v)))
        .collect::<Vec<_>>();
    let pf: Vec<Box<dyn LinearForm<BF>>> = vec![Box::new(Covector {
        vector: (0..sz as u64).map(BF::from).collect(),
    })];
    // Reset peak RIGHT before prove to get accurate per-run peak
    reset_peak();
    let t0 = Instant::now();
    let _ = p.prove(
        &mut ps,
        v.iter().map(|x| Cow::Borrowed(x.as_slice())).collect(),
        vec![Cow::Owned(w)],
        pf,
        Cow::Borrowed(vals.as_slice()),
    );
    let pm = t0.elapsed().as_secs_f64() * 1000.0;
    let pk = get_peak_mb();
    let proof = ps.proof();
    let t1 = Instant::now();
    let mut vs = VerifierState::new_std(&ds, &proof);
    let c = p.receive_commitment(&mut vs).unwrap();
    let wr = lfs
        .iter()
        .map(|w| w.as_ref() as &dyn LinearForm<BF>)
        .collect::<Vec<_>>();
    p.verify(&mut vs, &[&c], &vals).unwrap().verify(wr).unwrap();
    let vm = t1.elapsed().as_secs_f64() * 1000.0;
    (pm, pk, vm)
}

fn bench(sz: usize) -> (f64, f64, f64) {
    // Each bench: warmup then best of 3, with fresh peak each time
    let _ = run(sz);
    let mut b = run(sz);
    for _ in 0..2 {
        let r = run(sz);
        if r.0 < b.0 {
            b = r;
        }
    }
    b
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
    let d = 1usize << FF;
    let lo = 1usize << 12;
    let hi = 1usize << 18;

    // Generate ALL target sizes
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

    // For each target: compute smooth commit size
    println!("n,commit_size,prove_ms,peak_mb,verify_ms");
    for &n in &targets {
        let mut cs = next_smooth(n);
        while cs % d != 0 {
            cs = next_smooth(cs + 1);
        }
        let r = bench(cs);
        println!("{n},{cs},{:.3},{:.1},{:.3}", r.0, r.1, r.2);
        eprintln!(
            "  n={n:>7} -> {cs:>7}  prove={:.1}ms  mem={:.0}MB",
            r.0, r.1
        );
    }
}
