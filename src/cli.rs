use clap::{Parser, Subcommand};

use crate::commands::profile::ProfileArgs;

#[derive(Parser, Debug)]
#[command(name = "expanse", version, about = "Biobank-scale tandem repeat expansion screening")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Extract IRR (in-repeat-read) evidence and their mates from a CRAM/BAM
    /// into a small output CRAM/BAM for downstream repeat-expansion calling.
    Profile(ProfileArgs),
}
