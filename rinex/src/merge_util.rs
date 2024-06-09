//! RINEX File merging (combination)
use crate::prelude::Epoch;
use hifitime::EpochError;
use thiserror::Error;

/// Merge operation related error(s)
#[derive(Error, Debug)]
pub enum Error {
    #[error("cannot merge mixed absolute/relative phase antenna together")]
    AntexAbsoluteRelativeMismatch,
    #[error("cannot merge ionex based off different reference systems")]
    IonexReferenceMismatch,
    #[error("cannot merge ionex with different grid definition")]
    IonexMapGridMismatch,
    #[error("cannot merge ionex of different dimensions")]
    IonexMapDimensionsMismatch,
    #[error("cannot merge ionex where base radius differs")]
    IonexBaseRadiusMismatch,
    #[error("failed to retrieve system time for merge ops date")]
    HifitimeError(#[from] EpochError),
}

/*
 * Merges "TIME OF FIRST" special OBSERVATION header field
 */
pub(crate) fn merge_time_of_first_obs(lhs: &mut Option<Epoch>, rhs: &Option<Epoch>) {
    if lhs.is_none() {
        if let Some(rhs) = rhs {
            *lhs = Some(*rhs);
        }
    } else if let Some(rhs) = rhs {
        let tl = lhs.unwrap();
        *lhs = Some(std::cmp::min(tl, *rhs));
    }
}

/*
 * Merges "TIME OF LAST" special OBSERVATION header field
 */
pub(crate) fn merge_time_of_last_obs(lhs: &mut Option<Epoch>, rhs: &Option<Epoch>) {
    if lhs.is_none() {
        if let Some(rhs) = rhs {
            *lhs = Some(*rhs);
        }
    } else if let Some(rhs) = rhs {
        let tl = lhs.unwrap();
        *lhs = Some(std::cmp::max(tl, *rhs));
    }
}
