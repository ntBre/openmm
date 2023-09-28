#![feature(lazy_cell)]

use context::Context;
use integrators::Integrator;
use system::System;
use topology::Topology;

macro_rules! q {
    ($v:expr, $u:expr) => {
        Quantity::new($v, $u)
    };
}

mod element;

pub mod barostats;
pub mod context;
pub mod forcefield;
pub mod integrators;
pub mod pdb;
pub mod platform;
pub mod state;
pub mod system;
pub mod topology;

pub mod quantity;

/// Simulation provides a simplified API for running simulations with OpenMM and
/// reporting results.
///
/// A Simulation ties together various objects used for running a simulation: a
/// Topology, System, Integrator, and Context. To use it, you provide the
/// Topology, System, and Integrator, and it creates the Context automatically.
///
/// Simulation also maintains a list of "reporter" objects that record or
/// analyze data as the simulation runs, such as writing coordinates to files or
/// displaying structures on the screen. For example, the following line will
/// cause a file called "output.pdb" to be created, and a structure written to
/// it every 1000 time steps:
pub struct Simulation<I>
where
    I: Integrator,
{
    pub topology: Topology,
    pub system: System,
    pub integrator: I,
    pub context: Context,
}

impl<I> Simulation<I>
where
    I: Integrator,
{
    pub fn new(topology: Topology, system: System, mut integrator: I) -> Self {
        let context = Context::new(&system, &mut integrator);
        Self {
            topology,
            system,
            integrator,
            context,
        }
    }
}
