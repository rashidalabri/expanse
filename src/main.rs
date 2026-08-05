use clap::Parser;

use expanse::cli::{Cli, Commands};
use expanse::commands;

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();

    let result = match cli.command {
        Commands::Profile(args) => commands::profile::run(args),
    };

    if let Err(err) = result {
        log::error!("{err:#}");
        std::process::exit(1);
    }
}
