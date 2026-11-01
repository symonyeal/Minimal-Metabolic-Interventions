//! Minimal cut sets of metabolic networks: elementary modes by double description, minimal
//! hitting sets by MMCS, an independent checker and a census of search effort. File formats and
//! I/O live in the `mcs` binary.

pub mod bits;
pub mod cen;
pub mod chk;
pub mod dd;
pub mod hg;
pub mod net;
pub mod srch;
