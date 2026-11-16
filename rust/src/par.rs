//! par -- the parallel building blocks of the engine (rayon work stealing).
//!
//! Switches (environment): CCA2_THREADS=n sets the pool size (default: all hardware threads);
//! CCA2_THREADS=1 or CCA2_PAR=0 makes every helper run serially, in the original order, which is the
//! reference the parallel results are compared with.
//!
//! Contract of every helper: results come back IN INPUT ORDER and, when several items fail, the
//! error returned is the one the serial loop would have met first. The closures must be pure
//! (memo tables are caches, see memo.rs), so the output cannot depend on the schedule.
use crate::errors::R;
use rayon::prelude::*;
use std::sync::OnceLock;

/// Stack of a pool thread: the quadtree and arrangement recursions are deep (the main thread has 8 MB).
const STACK: usize = 128 << 20;

fn n_threads() -> usize {
    static N: OnceLock<usize> = OnceLock::new();
    *N.get_or_init(|| {
        if std::env::var("CCA2_PAR").map(|v| v == "0").unwrap_or(false) {
            return 1;
        }
        let n = std::env::var("CCA2_THREADS").ok().and_then(|v| v.parse::<usize>().ok()).filter(|&n| n >= 1);
        let n = n.unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1));
        if n > 1 {
            // ignore the error: a pool may already exist (tests, a host program); it is then used as is
            let _ = rayon::ThreadPoolBuilder::new().num_threads(n).stack_size(STACK).thread_name(|i| format!("cca2-{}", i)).build_global();
        }
        n
    })
}

/// True when the helpers may use more than one thread.
pub fn enabled() -> bool {
    n_threads() > 1
}

/// `items.iter().map(f)` collected in order, items in parallel. The first error in item order wins.
pub fn map<T: Sync, U: Send, F: Fn(&T) -> R<U> + Sync + Send>(items: &[T], f: F) -> R<Vec<U>> {
    if items.len() < 2 || !enabled() {
        return items.iter().map(f).collect();
    }
    let v: Vec<R<U>> = items.par_iter().map(f).collect();
    v.into_iter().collect()
}

/// Same for an index range 0..n.
pub fn map_range<U: Send, F: Fn(usize) -> R<U> + Sync + Send>(n: usize, f: F) -> R<Vec<U>> {
    if n < 2 || !enabled() {
        return (0..n).map(f).collect();
    }
    let v: Vec<R<U>> = (0..n).into_par_iter().map(f).collect();
    v.into_iter().collect()
}

/// True when reductions run as a balanced tree (CCA2_REDUCE=tree). The default is the left fold:
/// a maximum / minimum of QuaCons is the SAME FUNCTION whatever the association, but its face
/// decomposition is not canonical, so a tree gives an equal function with a different partition
/// (measured: conj fixture 35 has 71 faces folded, 67 as a tree). The default keeps the output
/// bit-identical to the serial engine.
pub fn tree_reduce() -> bool {
    static T: OnceLock<bool> = OnceLock::new();
    *T.get_or_init(|| std::env::var("CCA2_REDUCE").map(|v| v == "tree").unwrap_or(false))
}

/// Reduction of `items` (non-empty) with the associative `f`. Left fold by default (identical to the
/// serial engine). With CCA2_REDUCE=tree, a balanced tree: adjacent pairs are combined level by
/// level, so the depth is log2 n instead of n - 1 and the pairs of a level run in parallel.
pub fn reduce<T: Send + Sync, F: Fn(T, T) -> R<T> + Sync + Send>(mut items: Vec<T>, f: F) -> R<T> {
    assert!(!items.is_empty());
    if items.len() < 3 || !enabled() || !tree_reduce() {
        let mut it = items.into_iter();
        let mut g = it.next().unwrap();
        for x in it {
            g = f(g, x)?;
        }
        return Ok(g);
    }
    while items.len() > 1 {
        let odd = if items.len() % 2 == 1 { items.pop() } else { None };
        let mut pairs: Vec<(T, T)> = Vec::with_capacity(items.len() / 2);
        let mut it = items.into_iter();
        while let (Some(a), Some(b)) = (it.next(), it.next()) {
            pairs.push((a, b));
        }
        let res: Vec<R<T>> = pairs.into_par_iter().map(|(a, b)| f(a, b)).collect();
        items = res.into_iter().collect::<R<Vec<T>>>()?;
        if let Some(o) = odd {
            items.push(o);
        }
    }
    Ok(items.pop().unwrap())
}
