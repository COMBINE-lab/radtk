use clap::{Parser, Subcommand};

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
    let args = Cli::parse();
    println!("args = {:?}!", args);

    match args.command {
        Commands::Cat(cat_opts) => cat::cat(&cat_opts)?,
    }
    Ok(())
}
