use std::sync::{Arc, Mutex};
use crate::lookup_table_v2::LookupTable;
use crate::physics::{AnyPrimitive, LIGHT_SPEED};
use crate::traits::InitialModel;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

static UNIFORM_TEMPERATURE: f64 = 1e-6;

/// Jet propagating through a kilonova debris cloud and surrounding
/// relativistic envelop
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WindSelf {
    /// Rate of outflow of the wind
    pub wind_mass_outflow_rate: f64,

    /// Four velocity of wind
    pub wind_gamma_beta: f64,

    /// Rate of outflow of the flare
    #[serde(default)]
    pub flare_outflow_rate: f64,

    /// Four velocity of flare
    #[serde(default)]
    pub flare_gamma_beta: f64,

    /// Flare time
    #[serde(default)]
    pub flare_time: f64,

    /// Flare duration
    #[serde(default)]
    pub flare_duration: f64,

    #[serde(skip)]
    pub lookup_table: Arc<Mutex<Option<LookupTable<4>>>>,
}
/*
impl WindSelf {
    fn require_lookup_table(&self) {
        let mut self_table = self.lookup_table.as_ref().lock().unwrap();

        if self_table.is_none() {
            let filename = self.initial_data_table.as_ref().unwrap();
            let table = LookupTable::<4>::from_ascii_file(&filename).unwrap();
            *self_table = Some(table);
        }
    }
}
*/
// ============================================================================
impl InitialModel for WindSelf {
    fn validate(&self) -> anyhow::Result<()> {
        if self.wind_gamma_beta < 0.0 {
            anyhow::bail!("the wind four-velocity must be positive")
        }
        Ok(())
    }

    fn primitive_at(&self, coordinate: (f64, f64), t: f64) -> AnyPrimitive {
        // u: gamma-beta-c
        // v: beta-c
        // rho: comoving rest-mass density
        // Mdot = 4 pi r^2 rho u c

        if t >= self.flare_time && t <= self.flare_time + self.flare_duration {
            let r = coordinate.0;
            let u = if t >= self.flare_time && t < self.flare_time + 0.1 * self.flare_duration {
                self.wind_gamma_beta + (self.flare_gamma_beta - self.wind_gamma_beta) * (t - self.flare_time) / (0.1 * self.flare_duration)
            } else {
                self.wind_gamma_beta + ((t - self.flare_time - 0.1 * self.flare_duration) / (t - self.flare_time - self.flare_duration)).exp() * self.flare_gamma_beta
            };
            let mdot = if t >= self.flare_time && t <= self.flare_time + 0.9 * self.flare_duration {
                self.wind_mass_outflow_rate + self.flare_outflow_rate * (t - self.flare_time).powi(4) / (0.9 * self.flare_duration).powi(4)
            } else {
                self.flare_outflow_rate + (self.wind_mass_outflow_rate - self.flare_outflow_rate) * (t - self.flare_time - 0.9 * self.flare_duration) / (0.1 * self.flare_duration)
            };
            let rho = mdot / (4.0 * PI * r * r * u * LIGHT_SPEED);
            let p = rho * UNIFORM_TEMPERATURE;

            AnyPrimitive {
                velocity_r: u,
                velocity_q: 0.0,
                mass_density: rho,
                gas_pressure: p,
            }
        } else if t >= self.flare_time + 60.0 && t < self.flare_time + 60.0 + self.flare_duration {
            let r = coordinate.0;
            let u = if t >= self.flare_time + 60.0 && t < self.flare_time + 60.0 + 0.1 * self.flare_duration {
                self.wind_gamma_beta + (self.flare_gamma_beta - self.wind_gamma_beta) * (t - self.flare_time - 60.0) / (0.1 * self.flare_duration)
            } else {
                self.wind_gamma_beta + ((t - self.flare_time - 60.0 - 0.1 * self.flare_duration) / (t - self.flare_time - 60.0 - self.flare_duration)).exp() * self.flare_gamma_beta
            };
            let mdot = if t >= self.flare_time + 60.0 && t <= self.flare_time + 60.0 + 0.9 * self.flare_duration {
                self.wind_mass_outflow_rate + self.flare_outflow_rate * (t - self.flare_time - 60.0).powi(4) / (0.9 * self.flare_duration).powi(4)
            } else {
                self.flare_outflow_rate + (self.wind_mass_outflow_rate - self.flare_outflow_rate) * (t - self.flare_time - 60.0 - 0.9 * self.flare_duration) / (0.1 * self.flare_duration)
            };
            let rho = mdot / (4.0 * PI * r * r * u * LIGHT_SPEED);
            let p = rho * UNIFORM_TEMPERATURE;

            AnyPrimitive {
                velocity_r: u,
                velocity_q: 0.0,
                mass_density: rho,
                gas_pressure: p,
            }
        } else {
            let r = coordinate.0;
            let u = self.wind_gamma_beta;
            let rho = self.wind_mass_outflow_rate / (4.0 * PI * r * r * u * LIGHT_SPEED);
            let p = rho * UNIFORM_TEMPERATURE;
            AnyPrimitive {
                velocity_r: u,
                velocity_q: 0.0,
                mass_density: rho,
                gas_pressure: p,
            }
        }
    }

    fn scalar_at(&self, _coordinate: (f64, f64), _t: f64) -> f64 {
        0.0
    }
}
