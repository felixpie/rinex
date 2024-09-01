use crate::constants::Omega;
use dms_coordinates::DMS;

use anise::{
    astro::PhysicsResult,
    prelude::{Epoch, Frame, Orbit},
};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[cfg(feature = "qc")]
use maud::{html, Markup, Render};

#[derive(Copy, Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GroundPosition(Orbit);

impl GroundPosition {
    pub fn orbit(&self) -> Orbit {
        self.0
    }
    pub fn epoch(&self) -> Epoch {
        self.0.epoch
    }
    /// Builds [Self] from ECEF coordinates in km
    pub fn from_position_km(x_km: f64, y_km: f64, z_km: f64, t: Epoch, frame: Frame) -> Self {
        Self(Orbit::from_position(x_km, y_km, z_km, t, frame))
    }
    /// Builds [Self] from Geodetic angles in degrees
    pub fn from_geodetic(
        lat_deg: f64,
        long_deg: f64,
        h_km: f64,
        t: Epoch,
        frame: Frame,
    ) -> PhysicsResult<Self> {
        let angular_vel_deg_s = Omega::GPS_RAD_S.to_degrees();
        let orb = Orbit::try_latlongalt(lat_deg, long_deg, h_km, angular_vel_deg_s, t, frame)?;
        Ok(Self(orb))
    }
    /// Converts [Self] to ECEF coordinates in km
    pub fn to_position_km(&self) -> (f64, f64, f64) {
        let state = self.0.to_cartesian_pos_vel();
        (state[0], state[1], state[2])
    }
    /// Converts [Self] to geodetic angles in degrees
    pub fn to_geodetic(&self) -> PhysicsResult<(f64, f64, f64)> {
        self.0.latlongalt()
    }
    /// Returns position altitude
    pub fn altitude_km(&self) -> PhysicsResult<f64> {
        let (_, _, h_km) = self.to_geodetic()?;
        Ok(h_km)
    }
}

impl std::fmt::Display for GroundPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let (x_km, y_km, z_km) = self.to_position_km();
        write!(f, "x={}km, y={}km, z={}km", x_km, y_km, z_km)
    }
}

/*
 * RINEX compatible formatting
 */
impl std::fmt::UpperHex for GroundPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let (x_km, y_km, z_km) = self.to_position_km();
        let (x_m, y_m, z_m) = (x_km * 1.0E3, y_km * 1.0E3, z_km * 1.0E3);
        write!(f, "{:14.4}{:14.4}{:14.4}", x_m, y_m, z_m)
    }
}

#[cfg(feature = "qc")]
impl Render for GroundPosition {
    fn render(&self) -> Markup {
        let (x_km, y_km, z_km) = self.to_position_km();
        let (lat_deg, long_deg, h_km) = self.to_geodetic().unwrap_or_else(|_| (0.0, 0.0, 0.0));
        html! {
            table {
                tr {
                    th {
                        "ECEF"
                    }
                }
                tr {
                    th {
                        "X"
                    }
                    td {
                        (format!("{:.5} km", x_km))
                    }
                    th {
                        "Y"
                    }
                    td {
                        (format!("{:.5} km", y_km))
                    }
                    th {
                       "Z"
                    }
                    td {
                        (format!("{:.5} km", z_km))
                    }
                }
                tr {
                    th {
                        "GEO"
                    }
                }
                tr {
                    th {
                        "Latitude"
                    }
                    td {
                        (format!("{:.6}°", lat_deg))
                    }
                    th {
                        "Longitude"
                    }
                    td {
                        (format!("{:.6}°", long_deg))
                    }
                    th {
                        "Altitude"
                    }
                    td {
                        (format!("{:.5}m", h_km * 1.0E3))
                    }
                }
                tr {
                    th {
                        "DMS"
                    }
                    td {
                        (DMS::from_ddeg_latitude(lat_deg).to_string())
                    }
                    th {
                        "DMS"
                    }
                    td {
                        (DMS::from_ddeg_longitude(long_deg).to_string())
                    }
                }
            }
        }
    }
}
