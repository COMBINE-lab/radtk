use clap::{Parser, Subcommand};
use std::io::BufReader;

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

pub fn cat(cat_opts: &CatOpts) -> anyhow::Result<()> {
    use libradicl::rad_types;

    if let Some(fname) = cat_opts.inputs.first() {
        let f = std::fs::File::open(&fname)?;
        let mut ifile = BufReader::new(f);
        let p = libradicl::header::RadPrelude::from_bytes(&mut ifile)?;
        println!("{}", p.summary(Some(100))?);
    }
    Ok(())
}
