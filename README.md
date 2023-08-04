# openmm
safe Rust interface to my [openmm-sys](https://github.com/ntBre/openmm-sys)
package

[OpenMM](https://openmm.org/) is "a high-performance toolkit for molecular
simulation," written in C++. openmm-sys is a
[bindgen](https://github.com/rust-lang/rust-bindgen)-generated wrapper around
the C wrapper provided by OpenMM. This crate provides a high-level, safe Rust
interface around that C wrapper, modeled after the [OpenMM Python
API](http://docs.openmm.org/latest/api-python/).

## Limitations

Right now openmm-sys is hard-coded to use the path to OpenMM 8.0 on my machine.
If you're interested in using either of these crates, please let me know. It
should only take a one-line change to the [build
script](https://github.com/ntBre/openmm-sys/blob/7762614748d95b1e2f86b2e1729e0cbc456bdb36/build.rs#L6C53-L6C53),
but I don't want to pass around an environment variable to all of my
compilations while I'm the only one using it.

## Related projects

There is an existing
[openmm-sys](https://github.com/uob-positron-imaging-centre/openmm-rust) crate
on [crates.io](https://crates.io/crates/openmm-sys), but it does not seem to be
actively maintained, builds with OpenMM 7.4.2, and retrieves and builds OpenMM
automatically if you don't point it to your existing installation. That last
point is probably something my version should have eventually, but I'd rather
not do that for now.
