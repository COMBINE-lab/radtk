use clap::{Parser, Subcommand};
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};

mod cat;
use crate::cat::CatOpts;

/// testing out minimizer space suffix arrays
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(propagate_version = true)]
pub struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// concatenate the records in a sequence of RAD files
    Cat(CatOpts),
}

fn main() -> anyhow::Result<()> {
    // Check the `RUST_LOG` variable for the logger level and
    // respect the value found there. If this environment
    // variable is not set then set the logging level to
    // INFO.
    tracing_subscriber::registry()
        .with(fmt::layer())
        .with(
            EnvFilter::builder()
                .with_default_directive(LevelFilter::INFO.into())
                .from_env_lossy(),
        )
        .init();

    let args = Cli::parse();
    println!("args = {:?}!", args);

    match args.command {
        Commands::Cat(cat_opts) => cat::cat(&cat_opts)?,
    }
    Ok(())
}
