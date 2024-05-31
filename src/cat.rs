use anyhow::bail;
use clap::Parser;
use std::io::{BufReader, BufWriter, Write};
use tracing::{error, info, warn};

/// options relevant to building the minimizer space suffix array
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct CatOpts {
    /// ',' separated list of input RAD files
    #[arg(short, long, required = true, value_delimiter = ',')]
    inputs: Vec<std::path::PathBuf>,

    /// output RAD file
    #[arg(short, long, required = true)]
    output: std::path::PathBuf,
}

pub fn cat(cat_opts: &CatOpts) -> anyhow::Result<()> {
    if cat_opts.inputs.len() <= 1 {
        if let Some(input) = cat_opts.inputs.first() {
            warn!("You are attempting to concatenate a single input RAD file ({}) into a new output RAD file ({}); this operation does not make sense",
                input.display(),
                cat_opts.output.display()
            );
        } else {
            warn!("The input list of RAD files is empty; there's nothing to do.");
        }
        return Ok(());
    }

    let fname = cat_opts
        .inputs
        .first()
        .expect("input should contain multiple RAD files");
    let f = std::fs::File::open(&fname)?;
    let mut ifile = BufReader::new(f);
    let mut first_prelude = libradicl::header::RadPrelude::from_bytes(&mut ifile)?;
    let first_tag_map = first_prelude
        .file_tags
        .try_parse_tags_from_bytes(&mut ifile)?;
    drop(ifile);

    let mut total_chunks = first_prelude.hdr.num_chunks;

    for in_file in cat_opts.inputs.iter().skip(1) {
        let f = std::fs::File::open(&in_file)?;
        let mut ifile = BufReader::new(f);
        let new_prelude = libradicl::header::RadPrelude::from_bytes(&mut ifile)?;
        if new_prelude == first_prelude {
            total_chunks += new_prelude.hdr.num_chunks;
        } else {
            error!(
                "The prelude for ({}) is incompatible with the prelude for ({}); cannot proceed",
                cat_opts.inputs.first().unwrap().display(),
                in_file.display()
            );
            bail!("Incompatible input RAD files.");
        }
    }

    info!("All inputs had compatible preludes; continuing with merge!");
    info!("total chunks = {}", total_chunks);

    first_prelude.hdr.num_chunks = total_chunks;

    let ofile = std::fs::OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(&cat_opts.output)?;
    let mut owriter = BufWriter::new(ofile);

    // write the output prelude with the correct number of chunks.
    first_prelude
        .write(&mut owriter)
        .expect("cannot write output prelude to file");

    first_tag_map
        .write_values(&mut owriter)
        .expect("cannot write values of file-level tagl map to output file");

    for in_file in cat_opts.inputs.iter() {
        let f = std::fs::File::open(&in_file)?;
        let mut ifile = BufReader::new(f);
        let prelude = libradicl::header::RadPrelude::from_bytes(&mut ifile)?;
        let _tag_map = prelude.file_tags.try_parse_tags_from_bytes(&mut ifile)?;
        let copy_res = std::io::copy(&mut ifile, &mut owriter);
        if let Ok(copied_bytes) = copy_res {
            info!(
                "copied {} bytes of record chunks from {} into {}.",
                copied_bytes,
                in_file.display(),
                &cat_opts.output.display()
            );
        } else {
            bail!(
                "Failed to copy record chunks from {} to {}; error {:?}",
                in_file.display(),
                &cat_opts.output.display(),
                copy_res
            );
        }
    }

    Ok(())
}
