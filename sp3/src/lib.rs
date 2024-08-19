//! SP3 precise orbit file parser.
#![cfg_attr(docsrs, feature(doc_cfg))]
extern crate gnss_rs as gnss;

use itertools::Itertools;

#[cfg(feature = "qc")]
extern crate rinex_qc_traits as qc_traits;

use gnss::prelude::{Constellation, SV};
use hifitime::{Duration, Epoch, ParsingError as EpochParsingError, TimeScale};
use std::collections::BTreeMap;

use gnss_rs::constellation::ParsingError as ConstellationParsingError;
use std::str::FromStr;
use thiserror::Error;

#[cfg(feature = "processing")]
use qc_traits::processing::{
    Decimate, DecimationFilter, DecimationFilterType, FilterItem, MaskFilter, MaskOperand, Masking,
    Preprocessing,
};

use anise::prelude::{Almanac, Frame, Orbit};

#[cfg(test)]
mod tests;

mod header;
mod position;
mod reader;
mod velocity;
mod version;

#[cfg(docsrs)]
mod bibliography;

use header::{
    line1::{is_header_line1, Line1},
    line2::{is_header_line2, Line2},
};

use position::{position_entry, PositionEntry};
use velocity::{velocity_entry, VelocityEntry};

use reader::BufferedReader;
use std::io::BufRead;
use version::Version;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use std::path::Path;

pub mod prelude {
    pub use crate::{version::Version, DataType, Error, OrbitType, SP3};
    // Pub re-export
    pub use anise::prelude::{Almanac, Frame, Orbit};
    pub use gnss::prelude::{Constellation, SV};
    pub use hifitime::{Duration, Epoch, TimeScale};
}

fn file_descriptor(content: &str) -> bool {
    content.starts_with("%c")
}

fn sp3_comment(content: &str) -> bool {
    content.starts_with("/*")
}

fn end_of_file(content: &str) -> bool {
    content.eq("EOF")
}

fn new_epoch(content: &str) -> bool {
    content.starts_with("*  ")
}

#[derive(Default, Copy, Clone, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum DataType {
    #[default]
    Position,
    Velocity,
}

impl std::fmt::Display for DataType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Position => f.write_str("P"),
            Self::Velocity => f.write_str("V"),
        }
    }
}

impl std::str::FromStr for DataType {
    type Err = ParsingError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.eq("P") {
            Ok(Self::Position)
        } else if s.eq("V") {
            Ok(Self::Velocity)
        } else {
            Err(ParsingError::UnknownDataType(s.to_string()))
        }
    }
}

#[derive(Default, Copy, Clone, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum OrbitType {
    #[default]
    FIT,
    EXT,
    BCT,
    BHN,
    HLM,
}

impl std::fmt::Display for OrbitType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::FIT => f.write_str("FIT"),
            Self::EXT => f.write_str("EXT"),
            Self::BCT => f.write_str("BCT"),
            Self::BHN => f.write_str("BHN"),
            Self::HLM => f.write_str("HLM"),
        }
    }
}

impl std::str::FromStr for OrbitType {
    type Err = ParsingError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.eq("FIT") {
            Ok(Self::FIT)
        } else if s.eq("EXT") {
            Ok(Self::EXT)
        } else if s.eq("BCT") {
            Ok(Self::BCT)
        } else if s.eq("BHN") {
            Ok(Self::BHN)
        } else if s.eq("HLM") {
            Ok(Self::HLM)
        } else {
            Err(ParsingError::UnknownOrbitType(s.to_string()))
        }
    }
}

/*
 * Comments contained in file
 */
type Comments = Vec<String>;

/// [SP3Entry] indexer
#[derive(Debug, Clone, PartialEq, PartialOrd, Eq, Ord)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SP3Key {
    /// Spacecraft described as [SV]
    pub sv: SV,
    /// Epoch
    pub epoch: Epoch,
}

/// [SP3Entry] record file content, sorted per [SP3Key]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SP3Entry {
    /// Possible clock correction in microsecond with 1E-12 precision.
    /// Often omitted, files that only have orbital states are most common.
    pub clock: Option<f64>,
    /// Clock offset variation, with 0.1ns/s scaling with 1E-16 theoretical precision.
    /// Rarely present.
    pub clock_rate: Option<f64>,
    /// [Orbit] with state vector define in [km] with 1mm precision.
    /// Velocity might not be present, it depends on the input file
    pub orbit: Orbit,
}

impl SP3Entry {
    /// Builds new [SP3Entry] with ECEF position coordinates in km
    pub fn from_position(x_km: f64, y_km: f64, z_km: f64, t: Epoch, frame: Frame) -> Self {
        Self {
            clock: None,
            clock_rate: None,
            orbit: Orbit::from_position(x_km, y_km, z_km, t, frame),
        }
    }
    /// Builds new [SP3Entry] with ECEF position coordinates in km and velocity in km/s
    pub fn from_position_velocity(
        x_km: f64,
        y_km: f64,
        z_km: f64,
        vel_x_kms: f64,
        vel_y_kms: f64,
        vel_z_kms: f64,
    ) -> Self {
        let pos_vel = Vector6::new(x_km, y_km, z_km, vel_x_kms, vel_y_kms, vel_z_kms);
        Self {
            clock: None,
            clock_rate: None,
            orbit: Orbit::from_cartesian_pos_vel(pos_vel),
        }
    }
    /// Copies and returns [Self] with ECEF position coordinates in km
    pub fn with_position(&self, x_km: f64, y_km: f64, z_km: f64, t: Epoch, frame: Frame) -> Self {
        let mut s = self.clone();
        s.orbit = Orbit::from_position(x_km, y_km, z_km, t, frame);
        s
    }
    /// Copies and returns [Self] with velocity in km/s expressed in ECEF
    pub fn with_velocity(&self, x_km_s: f64, y_km_s: f64, z_km_s: f64) -> Self {
        let mut s = self.clone();
        let vel_km_s = Vector3::new(x_km_s, y_km_s, z_km_s);
        s.orbit.with_velocity_km_s(vel_km_s);
        s
    }
    /// Copies and returns [Self] with given clock offset
    pub fn with_clock_offset(&self, offset: f64) -> Self {
        let mut s = self.clone();
        s.clock = Some(offset);
        s
    }
    /// Copies and returns [Self] with given clock rate
    pub fn with_clock_rate(&self, rate: f64) -> Self {
        let mut s = self.clone();
        s.clock_rate = Some(rate);
        s
    }
}

#[derive(Default, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SP3 {
    /// File revision
    pub version: Version,
    /// Data Type used in this file.
    /// If DataType == Velocity, you know
    /// that velocities record will be provided.
    /// Otherwise, that is not garanteed and kind of rare.
    pub data_type: DataType,
    /// Coordinates system used in this file.
    pub coord_system: String,
    /// Type of Orbit contained in this file.
    pub orbit_type: OrbitType,
    /// Agency providing this data
    pub agency: String,
    /// Type of constellations encountered in this file.
    /// For example "GPS" means only GPS vehicles are present.
    pub constellation: Constellation,
    /// [TimeScale] that applies to all following [Epoch]
    pub time_scale: TimeScale,
    /// Initial week counter in [TimeScale]
    pub week_counter: (u32, f64),
    /// Initial MJD, in time_system
    pub mjd_start: (u32, f64),
    /// [`Epoch`]s where at least one position or clock offset is provided
    pub epoch: Vec<Epoch>,
    /// Returns sampling interval, ie., time between successive [`Epoch`]s.
    pub epoch_interval: Duration,
    /// Satellite Vehicles
    pub sv: Vec<SV>,
    /// File content are [SP3Entry]s sorted per [SP3Key]
    pub data: BTreeMap<SP3Key, SP3Entry>,
    /// File header comments, stored as is.
    pub comments: Comments,
}

#[derive(Debug, Error)]
pub enum Error {
    #[error("parsing error")]
    ParsingError(#[from] ParsingError),
    #[error("hifitime parsing error")]
    HifitimeParsingError(#[from] EpochParsingError),
    #[error("constellation parsing error")]
    ConstellationParsing(#[from] ConstellationParsingError),
    #[error("unknown or non supported revision \"{0}\"")]
    UnknownVersion(String),
    #[error("unknown data type \"{0}\"")]
    UnknownDataType(String),
    #[error("unknown orbit type \"{0}\"")]
    UnknownOrbitType(String),
    #[error("file i/o error")]
    DataParsingError(#[from] std::io::Error),
}

#[derive(Debug, Error)]
pub enum ParsingError {
    #[error("unknown or non supported revision \"{0}\"")]
    UnknownVersion(String),
    #[error("unknown data type \"{0}\"")]
    UnknownDataType(String),
    #[error("unknown orbit type \"{0}\"")]
    UnknownOrbitType(String),
    #[error("malformed header line #1")]
    MalformedH1,
    #[error("malformed header line #2")]
    MalformedH2,
    #[error("malformed %c line \"{0}\"")]
    MalformedDescriptor(String),
    #[error("failed to parse epoch year from \"{0}\"")]
    EpochYear(String),
    #[error("failed to parse epoch month from \"{0}\"")]
    EpochMonth(String),
    #[error("failed to parse epoch day from \"{0}\"")]
    EpochDay(String),
    #[error("failed to parse epoch hours from \"{0}\"")]
    EpochHours(String),
    #[error("failed to parse epoch minutes from \"{0}\"")]
    EpochMinutes(String),
    #[error("failed to parse epoch seconds from \"{0}\"")]
    EpochSeconds(String),
    #[error("failed to parse epoch milliseconds from \"{0}\"")]
    EpochMilliSeconds(String),
    #[error("failed to parse number of epochs \"{0}\"")]
    NumberEpoch(String),
    #[error("failed to parse week counter")]
    WeekCounter(String),
    #[error("failed to parse hifitime::Epoch")]
    Epoch,
    #[error("failed to parse sample rate from \"{0}\"")]
    EpochInterval(String),
    #[error("failed to parse mjd start \"{0}\"")]
    Mjd(String),
    #[error("failed to parse sv from \"{0}\"")]
    SV(String),
    #[error("failed to parse (x, y, or z) coordinates from \"{0}\"")]
    Coordinates(String),
    #[error("failed to parse clock data from \"{0}\"")]
    Clock(String),
}

/*
 * Parses hifitime::Epoch from standard format
 */
fn parse_epoch(content: &str, time_scale: TimeScale) -> Result<Epoch, ParsingError> {
    let y = u32::from_str(content[0..4].trim())
        .or(Err(ParsingError::EpochYear(content[0..4].to_string())))?;

    let m = u32::from_str(content[4..7].trim())
        .or(Err(ParsingError::EpochMonth(content[4..7].to_string())))?;

    let d = u32::from_str(content[7..10].trim())
        .or(Err(ParsingError::EpochDay(content[7..10].to_string())))?;

    let hh = u32::from_str(content[10..13].trim())
        .or(Err(ParsingError::EpochHours(content[10..13].to_string())))?;

    let mm = u32::from_str(content[13..16].trim())
        .or(Err(ParsingError::EpochMinutes(content[13..16].to_string())))?;

    let ss = u32::from_str(content[16..19].trim())
        .or(Err(ParsingError::EpochSeconds(content[16..19].to_string())))?;

    let _ss_fract = f64::from_str(content[20..27].trim()).or(Err(
        ParsingError::EpochMilliSeconds(content[20..27].to_string()),
    ))?;

    Epoch::from_str(&format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02} {}",
        y, m, d, hh, mm, ss, time_scale,
    ))
    .or(Err(ParsingError::Epoch))
}

impl SP3 {
    /// Parses given SP3 file, with possible seamless
    /// .gz decompression, if compiled with the "flate2" feature.
    pub fn from_path(path: &Path) -> Result<Self, Error> {
        let fullpath = path.to_string_lossy().to_string();
        Self::from_file(&fullpath)
    }
    /// See [Self::from_path]
    pub fn from_file(path: &str) -> Result<Self, Error> {
        let reader = BufferedReader::new(path)?;

        let mut version = Version::default();
        let mut data_type = DataType::default();

        let mut time_scale = TimeScale::default();
        let mut constellation = Constellation::default();
        let mut pc_count = 0_u8;

        let mut coord_system = String::from("Unknown");
        let mut orbit_type = OrbitType::default();
        let mut agency = String::from("Unknown");
        let mut week_counter = (0_u32, 0_f64);
        let mut epoch_interval = Duration::default();
        let mut mjd_start = (0_u32, 0_f64);

        let mut vehicles: Vec<SV> = Vec::new();
        let mut comments = Comments::new();
        let mut data = BTreeMap::<SP3Key, SP3Entry>::new();

        let mut epoch = Epoch::default();
        let mut epochs: Vec<Epoch> = Vec::new();

        for line in reader.lines() {
            let line = line.unwrap();
            let line = line.trim();
            if sp3_comment(line) {
                if line.len() > 4 {
                    comments.push(line[3..].to_string());
                }
                continue;
            }
            if end_of_file(line) {
                break;
            }
            if is_header_line1(line) && !is_header_line2(line) {
                let l1 = Line1::from_str(line)?;
                (version, data_type, coord_system, orbit_type, agency) = l1.to_parts();
            }
            if is_header_line2(line) {
                let l2 = Line2::from_str(line)?;
                (week_counter, epoch_interval, mjd_start) = l2.to_parts();
            }
            if file_descriptor(line) {
                if line.len() < 60 {
                    return Err(Error::ParsingError(ParsingError::MalformedDescriptor(
                        line.to_string(),
                    )));
                }

                if pc_count == 0 {
                    constellation = Constellation::from_str(line[3..5].trim())?;
                    time_scale = TimeScale::from_str(line[9..12].trim())?;
                }

                pc_count += 1;
            }
            if new_epoch(line) {
                epoch = parse_epoch(&line[3..], time_scale)?;
                epochs.push(epoch);
            }
            if position_entry(line) {
                if line.len() < 60 {
                    continue; // tolerates malformed positions
                }
                let entry = PositionEntry::from_str(line)?;
                let (sv, (x_km, y_km, z_km), clk) = entry.to_parts();

                //TODO : move this into %c config frame
                if !vehicles.contains(&sv) {
                    vehicles.push(sv);
                }
                // verify entry validity
                if x_km != 0.0_f64 && y_km != 0.0_f64 && z_km != 0.0_f64 {
                    let key = SP3Key { epoch, sv };
                    if let Some(e) = data.get_mut(&key) {
                        e.position = (x_km, y_km, z_km);
                    } else {
                        if let Some(clk) = clk {
                            data.insert(
                                key,
                                SP3Entry::from_position(x_km, y_km, z_km).with_clock_offset(clk),
                            );
                        } else {
                            data.insert(key, SP3Entry::from_position((x_km, y_km, z_km)));
                        }
                    }
                }
            }
            if velocity_entry(line) {
                if line.len() < 60 {
                    continue; // tolerates malformed velocities
                }
                let entry = VelocityEntry::from_str(line)?;
                let (sv, (vel_x, vel_y, vel_z), clk) = entry.to_parts();

                //TODO : move this into %c config frame
                if !vehicles.contains(&sv) {
                    vehicles.push(sv);
                }
                // verify entry validity
                if vel_x != 0.0_f64 && vel_y != 0.0_f64 && vel_z != 0.0_f64 {
                    let key = SP3Key { epoch, sv };
                    if let Some(e) = data.get_mut(&key) {
                        *e = e.with_velocity((vel_x, vel_y, vel_z));
                        if let Some(clk) = clk {
                            *e = e.with_clock_rate(clk);
                        }
                    } else {
                        if let Some(clk) = clk {
                            data.insert(
                                key,
                                SP3Entry::from_position(0.0, 0.0, 0.0, )
                                .with_clock_rate(clk),
                            );
                        } else {
                            data.insert(key, SP3Entry::from_position(0.0, 0.0, 0.0)));
                        }
                    }
                }
            }
        }
        Ok(Self {
            version,
            data_type,
            epoch: epochs,
            time_scale,
            constellation,
            coord_system,
            orbit_type,
            agency,
            week_counter,
            epoch_interval,
            mjd_start,
            sv: vehicles,
            data,
            comments,
        })
    }
    /// Returns a unique Epoch iterator where either
    /// Position or Clock data is provided.
    pub fn epoch(&self) -> impl Iterator<Item = Epoch> + '_ {
        self.epoch.iter().copied()
    }
    /// Returns total number of epoch
    pub fn nb_epochs(&self) -> usize {
        self.epoch.len()
    }
    /// Returns first epoch
    pub fn first_epoch(&self) -> Option<Epoch> {
        self.epoch.first().copied()
    }
    /// Returns last epoch
    pub fn last_epoch(&self) -> Option<Epoch> {
        self.epoch.last().copied()
    }
    /// Returns a unique [Constellation] iterator
    pub fn constellation(&self) -> impl Iterator<Item = Constellation> + '_ {
        self.sv().map(|sv| sv.constellation).unique()
    }
    /// Returns a unique [SV] iterator
    pub fn sv(&self) -> impl Iterator<Item = SV> + '_ {
        self.sv.iter().copied()
    }
    /// [Orbit] Iterator with state vectors expressed in km with 1mm precision.
    pub fn sv_orbit(&self) -> impl Iterator<Item = (Epoch, SV, Orbit)> + '_ {
        self.data.iter().map(|(k, v)| (k.epoch, k.sv, v.orbit))
    }
    /// Returns state vector Iterator, expressed in km with 1mm precision.
    pub fn sv_position_km(&self) -> impl Iterator<Item = (Epoch, SV, f64, f64, f64)> + '_ {
        self.sv_orbit().map(|(t, sv, orb)| {
            let cartesian = orb.to_cartesian_pos_vel();
            (t, sv, cartesian[0], cartesian[1], cartesian[2])
        })
    }
    /// Returns velocity vector Iterator, expressed in km
    pub fn sv_velocity_km_s(&self) -> impl Iterator<Item = (Epoch, SV, f64, f64, f64)> + '_ {
        self.sv_orbit().filter_map(|(t, sv, orb)| {
            if orb.vmag_km_s() > 0.0 {
                let cartesian = orb.to_cartesian_pos_vel();
                Some((t, sv, cartesian[3], cartesian[4], cartesian[5]))
            } else {
                None
            }
        })
    }
    /// [SV] attitude as (azimuth, elevation, range) Iterator.
    /// Units: degrees, degrees, km.
    pub fn sv_azimuth_elev_range(
        &self,
        rx_orbit: Orbit,
    ) -> impl Iterator<Item = (Epoch, SV, f64, f64, f64)> + '_ {
        self.sv_orbit().filter_map(|(t, sv, orb)| {
            if let Ok(azelrange) = azimuth_elevation_range(orb, rx_orbit) {
                let azim_deg = azelrange.azimuth_deg;
                let elev_deg = azelrange.elevation_deg;
                let range_km = azelrange.range;
                Some((t, sv, azim_deg, elev_deg, range_km))
            } else {
                None
            }
        })
    }
    /// Returns an Iterator over SV velocities estimates,
    /// in 10^-1 m/s with 0.1 um/s precision.
    pub fn sv_velocities(&self) -> impl Iterator<Item = (Epoch, SV, Vector3D)> + '_ {
        self.data.iter().filter_map(|(k, v)| {
            let velocity = v.velocity?;
            Some((k.epoch, k.sv, velocity))
        })
    }
    /// Returns an Iterator over Clock offsets with theoretical 1E-12 precision.
    pub fn sv_clock(&self) -> impl Iterator<Item = (Epoch, SV, f64)> + '_ {
        self.data.iter().filter_map(|(k, v)| {
            let clock = v.clock?;
            Some((k.epoch, k.sv, clock))
        })
    }
    /// Returns an Iterator over Clock offset variations, scaling is 0.1 ns/s and theoretical
    /// precision downto 0.1 fs/s precision.
    pub fn sv_clock_rate(&self) -> impl Iterator<Item = (Epoch, SV, f64)> + '_ {
        self.data.iter().filter_map(|(k, v)| {
            let rate = v.clock_rate?;
            Some((k.epoch, k.sv, rate))
        })
    }
    /// Returns an Iterator over [`Comments`] contained in this file
    pub fn comments(&self) -> impl Iterator<Item = &String> + '_ {
        self.comments.iter()
    }
}

#[cfg(feature = "qc")]
use qc_traits::{Merge, MergeError};

#[cfg(feature = "qc")]
impl Merge for SP3 {
    fn merge(&self, rhs: &Self) -> Result<Self, MergeError> {
        let mut s = self.clone();
        s.merge_mut(rhs)?;
        Ok(s)
    }
    fn merge_mut(&mut self, rhs: &Self) -> Result<(), MergeError> {
        if self.agency != rhs.agency {
            return Err(MergeError::DataProviderAgencyMismatch);
        }
        if self.time_scale != rhs.time_scale {
            return Err(MergeError::TimescaleMismatch);
        }
        if self.coord_system != rhs.coord_system {
            return Err(MergeError::ReferenceFrameMismatch);
        }
        if self.constellation != rhs.constellation {
            /*
             * Convert self to Mixed constellation
             */
            self.constellation = Constellation::Mixed;
        }
        // adjust revision
        if rhs.version > self.version {
            self.version = rhs.version;
        }
        // Adjust MJD start
        if rhs.mjd_start.0 < self.mjd_start.0 {
            self.mjd_start.0 = rhs.mjd_start.0;
        }
        if rhs.mjd_start.1 < self.mjd_start.1 {
            self.mjd_start.1 = rhs.mjd_start.1;
        }
        // Adjust week counter
        if rhs.week_counter.0 < self.week_counter.0 {
            self.week_counter.0 = rhs.week_counter.0;
        }
        if rhs.week_counter.1 < self.week_counter.1 {
            self.week_counter.1 = rhs.week_counter.1;
        }
        // update SV table
        for sv in &rhs.sv {
            if !self.sv.contains(sv) {
                self.sv.push(*sv);
            }
        }
        // update sampling interval (pessimistic)
        self.epoch_interval = std::cmp::max(self.epoch_interval, rhs.epoch_interval);
        // Merge new entries
        // and upgrade missing information (if possible)
        for (key, entry) in &rhs.data {
            if let Some(lhs_entry) = self.data.get_mut(key) {
                if let Some(clock) = entry.clock {
                    lhs_entry.clock = Some(clock);
                }
                if let Some(rate) = entry.clock_rate {
                    lhs_entry.clock_rate = Some(rate);
                }
                if lhs_entry.orbit.vmag_km_s() == 0.0 && entry.orbit.vmag_km_s() > 0.0 {
                    let pos_vel = entry.orbit.to_cartesian_pos_vel();
                    lhs_entry.orbit.with_velocity_km_s(Vector3::new(pos_vel[3], pos_vel[4], pos_vel[5]));
                }
            } else {
                if !self.epoch.contains(&key.epoch) {
                    self.epoch.push(key.epoch); // new epoch
                }
                self.data.insert(key.clone(), entry.clone()); // new entry
            }
        }
        self.epoch.sort(); // preserve chronological order
        Ok(())
    }
}

#[cfg(feature = "processing")]
impl Preprocessing for SP3 {}

#[cfg(feature = "processing")]
impl Masking for SP3 {
    fn mask(&self, f: &MaskFilter) -> Self {
        let mut s = self.clone();
        s.mask_mut(&f);
        s
    }
    fn mask_mut(&mut self, f: &MaskFilter) {
        match f.operand {
            MaskOperand::Equals => match &f.item {
                FilterItem::EpochItem(epoch) => {
                    self.data.retain(|k, _| k.epoch == *epoch);
                },
                FilterItem::SvItem(svs) => {
                    self.data.retain(|k, _| svs.contains(&k.sv));
                },
                FilterItem::ConstellationItem(constells) => {
                    let mut broad_sbas_filter = false;
                    for c in constells {
                        broad_sbas_filter |= *c == Constellation::SBAS;
                    }
                    self.data.retain(|k, _| {
                        if broad_sbas_filter {
                            k.sv.constellation.is_sbas() || constells.contains(&k.sv.constellation)
                        } else {
                            constells.contains(&k.sv.constellation)
                        }
                    });
                },
                _ => {}, // does not apply
            },
            MaskOperand::NotEquals => match &f.item {
                FilterItem::EpochItem(epoch) => {
                    self.data.retain(|k, _| k.epoch != *epoch);
                },
                FilterItem::SvItem(svs) => {
                    self.data.retain(|k, _| !svs.contains(&k.sv));
                },
                FilterItem::ConstellationItem(constells) => {
                    self.data
                        .retain(|k, _| !constells.contains(&k.sv.constellation));
                },
                _ => {}, // does not apply
            },
            MaskOperand::GreaterThan => match &f.item {
                FilterItem::EpochItem(epoch) => {
                    self.data.retain(|k, _| k.epoch > *epoch);
                },
                FilterItem::SvItem(svs) => {
                    self.data.retain(|k, _| {
                        let mut retain = false;
                        for sv in svs {
                            if k.sv.constellation == sv.constellation {
                                retain = k.sv.prn > sv.prn
                            } else {
                                retain = false
                            }
                        }
                        retain
                    });
                },
                _ => {}, // does not apply
            },
            MaskOperand::GreaterEquals => match &f.item {
                FilterItem::EpochItem(epoch) => {
                    self.data.retain(|k, _| k.epoch >= *epoch);
                },
                FilterItem::SvItem(svs) => {
                    self.data.retain(|k, _| {
                        let mut retain = false;
                        for sv in svs {
                            if k.sv.constellation == sv.constellation {
                                retain = k.sv.prn >= sv.prn
                            } else {
                                retain = false
                            }
                        }
                        retain
                    });
                },
                _ => {}, // does not apply
            },
            MaskOperand::LowerThan => match &f.item {
                FilterItem::EpochItem(epoch) => {
                    self.data.retain(|k, _| k.epoch < *epoch);
                },
                FilterItem::SvItem(svs) => {
                    self.data.retain(|k, _| {
                        let mut retain = false;
                        for sv in svs {
                            if k.sv.constellation == sv.constellation {
                                retain = k.sv.prn < sv.prn
                            } else {
                                retain = false
                            }
                        }
                        retain
                    });
                },
                _ => {}, // does not apply
            },
            MaskOperand::LowerEquals => match &f.item {
                FilterItem::EpochItem(epoch) => {
                    self.data.retain(|k, _| k.epoch <= *epoch);
                },
                FilterItem::SvItem(svs) => {
                    self.data.retain(|k, _| {
                        let mut retain = false;
                        for sv in svs {
                            if k.sv.constellation == sv.constellation {
                                retain = k.sv.prn <= sv.prn
                            } else {
                                retain = false
                            }
                        }
                        retain
                    });
                },
                _ => {}, // does not apply
            },
        }
    }
}

#[cfg(feature = "processing")]
impl Decimate for SP3 {
    fn decimate(&self, f: &DecimationFilter) -> Self {
        let mut s = self.clone();
        s.decimate_mut(&f);
        s
    }
    fn decimate_mut(&mut self, f: &DecimationFilter) {
        if f.item.is_some() {
            todo!("targetted decimation not supported yet");
        }
        match f.filter {
            DecimationFilterType::Modulo(r) => {
                self.epoch_interval = self.epoch_interval * r as f64;
                let mut i = 0;
                self.data.retain(|_, _| {
                    let retained = (i % r) == 0;
                    i += 1;
                    retained
                });
            },
            DecimationFilterType::Duration(interval) => {
                self.epoch_interval = interval;
                let mut last_retained = Option::<Epoch>::None;
                self.data.retain(|k, _| {
                    if let Some(last) = last_retained {
                        let dt = k.epoch - last;
                        if dt >= interval {
                            last_retained = Some(k.epoch);
                            true
                        } else {
                            false
                        }
                    } else {
                        last_retained = Some(k.epoch);
                        true
                    }
                });
            },
        }
    }
}
