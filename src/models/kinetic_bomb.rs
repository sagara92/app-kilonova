use std::f64::consts::PI;
use serde::{Serialize, Deserialize};
use crate::traits::InitialModel;
use crate::physics::{AnyPrimitive, LIGHT_SPEED};

const UNIFORM_TEMPERATURE: f64 = 1e-3;




/**
 * Explosion in a horizontally stratified external medium
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct KineticBomb {
    pub external_medium_density: f64,
    pub launch_radius: f64,
    pub shell_thickness: f64,
    pub kinetic_energy: f64,
    pub shell_mass: f64,
}




// ============================================================================
impl KineticBomb {


    /**
     * Return the radial extent (in cm) of the shell at time t.
     */
    fn shell_extent(&self, t: f64) -> std::ops::Range<f64> {
        let r_outer_shell_surface = self.launch_radius + self.shell_velocity() * t;
        let r_inner_shell_surface = self.launch_radius + self.shell_velocity() * (t - self.shell_duration());
        r_inner_shell_surface..r_outer_shell_surface
    }


    /**
     * The velocity (in cm/s) the shell moves at (computed from the shell mass
     * and kinetic energy). Note this is expression assumes the shell is
     * sub-relativistic.
     */
    fn shell_velocity(&self) -> f64 {
        (2.0 * self.kinetic_energy / self.shell_mass).sqrt()
    }


    /**
     * The duration (in s) during which the shell is emerging from the inner
     * boundary.
     */
    fn shell_duration(&self) -> f64 {
        self.shell_thickness / self.shell_velocity()
    }
}




// ============================================================================
impl InitialModel for KineticBomb {

    fn validate(&self) -> anyhow::Result<()> {
        if self.shell_velocity() > 0.25 * LIGHT_SPEED {
            anyhow::bail!{"
             The shell is moving faster (v/c = {}) than 0.25 c, but
             this problem assumes Newtonian expressions for the
             kinetic energy. Consider reducing the kinetic energy or
             increasing the shell mass.", self.shell_velocity() / LIGHT_SPEED}
        } else {
            Ok(())
        }
    }

    fn primitive_at(&self, coordinate: (f64, f64), t: f64) -> AnyPrimitive {
        let (r, _q) = coordinate;

        if self.shell_extent(t).contains(&r) {
            let mdot = self.shell_mass / self.shell_duration();
            let v = self.shell_velocity();
            let d = mdot / (4.0 * PI * r * r * v);
            let p = d * UNIFORM_TEMPERATURE;

            AnyPrimitive {
                velocity_r: v / LIGHT_SPEED,
                velocity_q: 0.0,
                mass_density: d,
                gas_pressure: p,
            }            
        } else {
            let d0 = self.external_medium_density;
            let d = d0 * (r / self.launch_radius).powi(2);
            let p = d * UNIFORM_TEMPERATURE;

            AnyPrimitive {
                velocity_r: 0.0,
                velocity_q: 0.0,
                mass_density: d,
                gas_pressure: p,
            }            
        }
    }

    fn scalar_at(&self, coordinate: (f64, f64), t: f64) -> f64 {
        let (r, _q) = coordinate;
        if self.shell_extent(t).contains(&r) {
            1.0
        } else {
            0.0
        }
    }
}
