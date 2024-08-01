use clap::Parser;
use scroll::Pread;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use tracing::info;

/// options relevant to building the minimizer space suffix array
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct SplitOpts {
    /// input RAD file to split
    #[arg(short, long, required = true)]
    input: std::path::PathBuf,

    /// approximate number of reads in each sub-RAD file
    /// (Note: This is approximate because file chunks will not be split
    /// and input chunks will reside entirely within a single output file).
    #[arg(short, long, required = true)]
    num_reads: usize,

    /// output prefix
    #[arg(short, long, required = true)]
    output_prefix: std::path::PathBuf,

    /// be quiet (no progress bar or standard output messages)
    #[arg(short, long)]
    quiet: bool,
}

// TODO: There should be a "chunk-type-agnostic" read header function in `libradicl`
// add this.
fn read_chunk_header<F: std::io::BufRead>(f: &mut F) -> anyhow::Result<(u32, u32)> {
    let mut buf = [0u8; 8];
    f.read_exact(&mut buf)?;
    let nbytes = buf.pread::<u32>(0)?;
    let nrec = buf.pread::<u32>(4)?;
    Ok((nbytes, nrec))
}

fn process_file<F: std::io::BufRead + std::io::Seek>(
    f: &mut F,
    total_size: u64,
    in_prelude: &mut libradicl::header::RadPrelude,
    split_opts: &SplitOpts,
) -> anyhow::Result<()> {
    let mut file_ctr = 0_usize;
    let mut rec_in_current_output = 0_usize;
    let tag_map = in_prelude.file_tags.try_parse_tags_from_bytes(f)?;
    in_prelude.hdr.num_chunks = 0;

    let out_name_base = split_opts.output_prefix.clone();
    let mut out_name = out_name_base.clone();
    out_name.set_extension(format!("{}.rad", file_ctr));
    if out_name.exists() {
        std::fs::remove_file(&out_name)?;
    }

    let mut out_writer = BufWriter::new(std::fs::File::create(out_name.clone())?);
    let mut chunk_buf = Vec::<u8>::new();

    let current_offset = f.stream_position().expect("should be able to seek");
    let remaining = total_size.saturating_sub(current_offset);

    let pbar = indicatif::ProgressBar::new(remaining);
    pbar.set_draw_target(indicatif::ProgressDrawTarget::stderr_with_hz(5));
    if split_opts.quiet {
        pbar.set_draw_target(indicatif::ProgressDrawTarget::hidden());
    } else {
        pbar.set_style(
            indicatif::ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes}",
            )
                .unwrap()
                .progress_chars("#>-"),
        );
    }
    // write the header
    in_prelude.write(&mut out_writer)?;
    tag_map.write_values(&mut out_writer)?;

    while libradicl::utils::has_data_left(f).expect("encountered error reading input file") {
        let (num_bytes, num_rec) = read_chunk_header(f)?;

        let num_new_rec = num_rec as usize;
        if rec_in_current_output > 0
            && (rec_in_current_output + num_new_rec >= split_opts.num_reads)
        {
            // finish writing the old file.
            out_writer.flush()?;

            // create the new file
            file_ctr += 1;
            out_name.clone_from(&out_name_base);
            out_name.set_extension(format!("{}.rad", file_ctr));
            if out_name.exists() {
                std::fs::remove_file(&out_name)?;
            }
            out_writer = BufWriter::new(std::fs::File::create(out_name.clone())?);

            // write the header
            in_prelude.write(&mut out_writer)?;
            tag_map.write_values(&mut out_writer)?;

            // reset rec counter
            rec_in_current_output = 0;
        }
        rec_in_current_output += num_new_rec;
        // copy the chunk
        // first write the header
        out_writer.write_all(&num_bytes.to_le_bytes())?;
        out_writer.write_all(&num_rec.to_le_bytes())?;
        // copy the rest of the chunk
        chunk_buf.resize((num_bytes - 8) as usize, 0);
        f.read_exact(chunk_buf.as_mut_slice())?;
        std::io::copy(&mut &chunk_buf[..], &mut out_writer)?;
        pbar.inc(num_bytes as u64);
    }
    out_writer.flush()?;
    pbar.finish();
    if !split_opts.quiet {
        info!("generated {} output RAD files", file_ctr + 1);
    }
    Ok(())
}

pub fn split(split_opts: &SplitOpts) -> anyhow::Result<()> {
    let fname = split_opts.input.clone();

    let md = std::fs::metadata(&fname)?;
    let f = std::fs::File::open(fname)?;
    let file_size = md.len();
    let mut ifile = BufReader::new(f);
    let mut in_prelude = libradicl::header::RadPrelude::from_bytes(&mut ifile)?;
    process_file(&mut ifile, file_size, &mut in_prelude, split_opts)
}
