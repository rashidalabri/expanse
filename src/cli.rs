use clap::{Parser, Subcommand};

use crate::commands::merge::MergeArgs;
use crate::commands::profile::ProfileArgs;
use crate::commands::sinks::SinksArgs;

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

    /// Merge `profile --summary` outputs from many samples (given by a
    /// manifest) into one combined summary, with per-sample IRR counts.
    Merge(MergeArgs),

    /// Scan an entire CRAM/BAM for in-repeat reads (IRRs) and report where
    /// they cluster as a BED file, for use as a `profile --sink-bed` /
    /// `--exclude-bed` input.
    Sinks(SinksArgs),
}
