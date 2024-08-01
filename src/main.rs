use clap::{Parser, Subcommand};
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};

mod cat;
use crate::cat::CatOpts;
mod view;
use crate::view::ViewOpts;
mod split;
use crate::split::SplitOpts;

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
    /// view a text representation of the RAD file
    View(ViewOpts),
    /// split an input RAD file into multiple output files
    Split(SplitOpts),
}

fn main() -> anyhow::Result<()> {
    // Check the `RUST_LOG` variable for the logger level and
    // respect the value found there. If this environment
    // variable is not set then set the logging level to
    // INFO.
    tracing_subscriber::registry()
        .with(fmt::layer().with_writer(std::io::stderr))
        .with(
            EnvFilter::builder()
                .with_default_directive(LevelFilter::INFO.into())
                .from_env_lossy(),
        )
        .init();

    let args = Cli::parse();

    match args.command {
        Commands::Cat(cat_opts) => cat::cat(&cat_opts)?,
        Commands::View(view_opts) => view::view(&view_opts)?,
        Commands::Split(split_opts) => split::split(&split_opts)?,
    }
    Ok(())
}
