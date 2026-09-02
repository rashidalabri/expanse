use clap::Parser;

use expanse::cli::{Cli, Commands};
use expanse::commands;

fn main() {
    let cli = Cli::parse();

    let default_level = match cli.verbose {
        0 => "info",
        1 => "debug",
        _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(default_level))
        .init();

    let result = match cli.command {
        Commands::Profile(args) => commands::profile::run(args),
        Commands::Merge(args) => commands::merge::run(args),
    };

    if let Err(err) = result {
        log::error!("{err:#}");
        std::process::exit(1);
    }
}
