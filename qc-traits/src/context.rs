//! Traits required to form GNSS processing contexts
use thiserror::Error;

/// Error during [Merge]ing process
#[derive(Error, Debug)]
pub enum MergeError {
    #[error("non supported file format")]
    NonSupportedFileType,
    #[error("type mismatch: can't merge B into A")]
    FileTypeMismatch,
    #[error("merge error: {0}")]
    Custom(String),
}

/// The [Merge] Trait is describes context extension,
/// by appending data into the context. See the
/// implementation in RINEX or SP3 for examples.
pub trait Merge {
    /// Merge "rhs" dataset into self, to form an extended dataset.
    fn merge(&self, rhs: &Self) -> Result<Self, MergeError>
    where
        Self: Sized;
    /// Merge with mutable access, avoiding one memcopy.
    fn merge_mut(&mut self, rhs: &Self) -> Result<(), MergeError>;
}

/// Utilities that make [Merge]ing process easier
pub mod util {
    use std::hash::Hash;
    use std::collections::HashMap;
    /// Merges Option<B> into Option<A> only if A is None.
    pub fn merge_mut_option<T: Clone>(lhs: &mut Option<T>, rhs: &Option<T>) {
        if lhs.is_none() {
            if let Some(rhs) = rhs {
                *lhs = Some(rhs.clone());
            }
        }
    }
    /// Append one vec<B> into Vec<A> whathever the values
    pub fn merge_mut_vec<T: Clone>(lhs: &mut Vec<T>, rhs: &Vec<T>) {
        for item in rhs {
            lhs.push(item.clone());
        }
    }
    /// Merge B into A maintaining a unique list.
    pub fn merge_mut_unique_vec<T: Clone + PartialEq>(lhs: &mut Vec<T>, rhs: &Vec<T>) {
        for item in rhs {
            if !lhs.contains(item) {
                lhs.push(item.clone());
            }
        }
    }
    /// Merge 2D Map B into A maintaining unique keys
    pub fn merge_mut_unique_map2d<K: PartialEq + Eq + Hash + Clone, V: Clone + PartialEq>(
        lhs: &mut HashMap<K, Vec<V>>,
        rhs: &HashMap<K, Vec<V>>,
    ) {
        for (k, values) in rhs.iter() {
            if let Some(vvalues) = lhs.get_mut(k) {
                for value in values {
                    if !vvalues.contains(value) {
                        vvalues.push(value.clone());
                    }
                }
            } else {
                lhs.insert(k.clone(), values.clone());
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::util::*;
    #[test]
    fn test_merge_mut_option() {
        let mut a: Option::<u32> = Some(0);
        let b: Option::<u32> = Some(1);
        merge_mut_option(&mut a, &b);
        assert_eq!(a, Some(0));

        let mut a = Option::<u32>::None;
        merge_mut_option(&mut a, &b);
        assert_eq!(a, Some(1));
    }
}
