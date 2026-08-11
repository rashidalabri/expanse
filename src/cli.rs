use clap::{Parser, Subcommand};

use crate::commands::find_sinks::FindSinksArgs;
use crate::commands::profile::ProfileArgs;

#[derive(Parser, Debug)]
#[command(name = "expanse", version, about = "Biobank-scale tandem repeat expansion screening")]
pub struct Cli {
    /// Increase logging verbosity (-v for debug, -vv for trace).
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, global = true)]
    pub verbose: u8,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Extract IRR (in-repeat-read) evidence and their mates from a CRAM/BAM
    /// into a small output CRAM/BAM for downstream repeat-expansion calling.
    Profile(ProfileArgs),

    /// Scan an entire CRAM/BAM/SAM front-to-back (no index required) for
    /// in-repeat reads (IRRs), and write a BED of the "sink" regions they
    /// cluster into.
    FindSinks(FindSinksArgs),
}
