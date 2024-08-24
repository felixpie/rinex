//! const value of Gnss

use crate::prelude::{Constellation, SV};

pub(crate) struct GM;

impl GM {
    pub const GPS: f64 = 3.9860050E14;
    pub const BDS: f64 = 3.986004418E14;
    pub const GLO: f64 = 3.9860044E14;
    pub const GAL: f64 = 3.986004418E14;
}

/// Angular velocity in rad/s
pub(crate) struct Omega;

impl Omega {
    pub const GPS_RAD_S: f64 = 7.2921151467E-5;
    pub const BDS_RAD_S: f64 = 7.292115E-5;
    pub const GLO_RAD_S: f64 = 7.292115E-5;
    pub const GAL_RAD_S: f64 = 7.2921151467E-5;
}

/// - 2 * sqrt(gm) / c / c
#[allow(dead_code)]
pub(crate) struct DtrF;
#[allow(dead_code)]
impl DtrF {
    pub const GPS: f64 = -0.000000000444280763339306;
    pub const BDS: f64 = -0.00000000044428073090439775;
    pub const GAL: f64 = -0.00000000044428073090439775;
}

pub(crate) struct MaxIterNumber;

#[allow(dead_code)]
impl MaxIterNumber {
    /// Maximum number of iterations to calculate the anastomosis angle
    pub const KEPLER: u8 = 30;
}

/// Const values used in RINEX crate
pub(crate) struct Constants;

impl Constants {
    // earth
    pub const fn gm(sv: SV) -> f64 {
        match sv.constellation {
            Constellation::BeiDou => GM::BDS,
            Constellation::Galileo => GM::GAL,
            Constellation::Glonass => GM::GLO,
            _ => GM::GPS,
        }
    }
    /// Earth rotation rate
    pub const fn omega_rad_s(sv: SV) -> f64 {
        match sv.constellation {
            Constellation::BeiDou => Omega::BDS_RAD_S,
            Constellation::Galileo => Omega::GAL_RAD_S,
            Constellation::Glonass => Omega::GLO_RAD_S,
            _ => Omega::GPS_RAD_S,
        }
    }
    ///  Auxiliary Quantities for Calculating Relativistic Effects in Clock Correction
    pub const fn dtr_f(sv: SV) -> f64 {
        match sv.constellation {
            Constellation::GPS => DtrF::GPS,
            Constellation::BeiDou => DtrF::BDS,
            Constellation::Galileo => DtrF::GAL,
            _ => DtrF::GPS,
        }
    }
}
