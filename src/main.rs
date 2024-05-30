use clap::{Parser, Subcommand};

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

/// options relevant to building the minimizer space suffix array
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct CatOpts {
    /// ',' separated list of input RAD files
    #[arg(short, long, value_delimiter = ',')]
    inputs: Vec<std::path::PathBuf>,

    /// output RAD file
    #[arg(short, long)]
    output: std::path::PathBuf,
}

fn main() {
    let args = Cli::parse();
    println!("args = {:?}!", args);
}
