//! memo -- the shared memo tables. Every table is a PURE cache (the value is a function of the key,
//! and a size cap that empties it already existed), so sharing one table between threads can change
//! the running time but never a result. Sharded `RwLock<HashMap>`: lookups take a read lock on one
//! shard, inserts a write lock on one shard; the value is computed OUTSIDE any lock, so two threads
//! may compute the same entry once each (the second insert overwrites with an equal value).
use std::borrow::Borrow;
use std::collections::hash_map::RandomState;
use std::collections::HashMap;
use std::hash::{BuildHasher, Hash};
use std::sync::RwLock;

const SHARDS: usize = 32;

pub struct Memo<K, V> {
    hasher: RandomState,
    shard_cap: usize,
    shards: Vec<RwLock<HashMap<K, V>>>,
}

impl<K: Hash + Eq, V: Clone> Memo<K, V> {
    /// `cap`: the table is emptied (per shard, cap / SHARDS entries) when it grows past this.
    pub fn new(cap: usize) -> Self {
        Memo { hasher: RandomState::new(), shard_cap: (cap / SHARDS).max(1), shards: (0..SHARDS).map(|_| RwLock::new(HashMap::new())).collect() }
    }

    fn shard<Q: ?Sized + Hash>(&self, k: &Q) -> &RwLock<HashMap<K, V>> {
        &self.shards[(self.hasher.hash_one(k) as usize) % SHARDS]
    }

    pub fn get<Q: ?Sized + Hash + Eq>(&self, k: &Q) -> Option<V>
    where
        K: Borrow<Q>,
    {
        self.shard(k).read().unwrap().get(k).cloned()
    }

    pub fn insert(&self, k: K, v: V) {
        let mut m = self.shard(&k).write().unwrap();
        if m.len() >= self.shard_cap {
            m.clear();
        }
        m.insert(k, v);
    }

    pub fn clear(&self) {
        for s in &self.shards {
            s.write().unwrap().clear();
        }
    }
}
