use anyhow::bail;
use clap::{Parser, ValueEnum};
use libradicl::record::{PiscemBulkReadRecord, PiscemBulkRecordContext};
use std::fs::File;
use std::io;
use std::io::{BufReader, Write};
use tracing::{error, info, warn};

#[derive(Clone, Debug, PartialEq, ValueEnum)]
pub enum RadFileType {
    Bulk,
    SingleCell,
    Unknown,
}

/// options relevant to building the minimizer space suffix array
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct ViewOpts {
    /// ',' separated list of input RAD files
    #[arg(short, long, required = true)]
    input: std::path::PathBuf,

    /// output RAD file
    #[arg(short, long)]
    output: Option<std::path::PathBuf>,

    /// the type of input RAD file
    #[arg(short, long)]
    rad_type: RadFileType,

    /// skip printing the header and file-level tags (i.e. only print the mapping reacords)
    #[arg(long)]
    no_header: bool,
}

pub fn view(view_opts: &ViewOpts) -> anyhow::Result<()> {
    if view_opts.rad_type == RadFileType::Unknown {
        error!("Unknown file type not yet supported");
        bail!("Unknown file type not yet supported");
    }

    let mut output_stream: Box<dyn Write> = match view_opts.output {
        Some(ref path) => File::open(path).map(|f| Box::new(f) as Box<dyn Write>)?,
        None => Box::new(io::stdout()),
    };

    let f = std::fs::File::open(&view_opts.input)?;
    let mut ifile = BufReader::new(f);
    let prelude = libradicl::header::RadPrelude::from_bytes(&mut ifile)?;
    let _tag_map = prelude.file_tags.try_parse_tags_from_bytes(&mut ifile)?;

    writeln!(output_stream, "{{")?;
    if !view_opts.no_header {
        writeln!(output_stream, "\"rad_header\" : {{")?;
        writeln!(output_stream, "\"is_paired\" : {},", prelude.hdr.is_paired)?;
        writeln!(output_stream, "\"ref_count\" : {},", prelude.hdr.ref_count)?;

        writeln!(output_stream, "\"refs\" : [")?;
        for (i, rn) in prelude.hdr.ref_names.iter().enumerate() {
            write!(output_stream, "\"{}\"", rn)?;
            if i < (prelude.hdr.ref_count - 1) as usize {
                writeln!(output_stream, ",")?;
            }
        }
        writeln!(output_stream, "],")?;

        writeln!(output_stream, "\"num_chunks\" : {}", prelude.hdr.num_chunks)?;
        writeln!(output_stream, "}},")?;

        writeln!(output_stream, "\"tag_descriptions\" : {{")?;

        writeln!(output_stream, "\t\"file_tag_desc\" : {{")?;
        writeln!(
            output_stream,
            "\t\t\"label\" : \"{:?}\",",
            prelude.file_tags.label
        )?;

        writeln!(output_stream, "\t\t\"tag_desc\" : [")?;
        for (i, td) in prelude.file_tags.tags.iter().enumerate() {
            writeln!(
                output_stream,
                "\t\t\t{{\n\t\t\t \"name\" : \"{}\",",
                td.name
            )?;
            write!(
                output_stream,
                "\t\t\t \"desc\" : \"{:?}\"\n\t\t\t}}",
                td.typeid
            )?;
            if i < (prelude.file_tags.tags.len() - 1) as usize {
                writeln!(output_stream, ",")?;
            } else {
                write!(output_stream, "\n")?;
            }
        }
        writeln!(output_stream, "\t\t]")?;
        writeln!(output_stream, "\t}},")?;

        writeln!(output_stream, "\t\"read_tag_desc\" : {{")?;
        writeln!(
            output_stream,
            "\t\t\"label\" : \"{:?}\",",
            prelude.read_tags.label
        )?;

        writeln!(output_stream, "\t\t\"tag_desc\" : [")?;
        for (i, td) in prelude.read_tags.tags.iter().enumerate() {
            writeln!(
                output_stream,
                "\t\t\t{{\n\t\t\t \"name\" : \"{}\",",
                td.name
            )?;
            write!(
                output_stream,
                "\t\t\t \"desc\" : \"{:?}\"\n\t\t\t}}",
                td.typeid
            )?;
            if i < (prelude.read_tags.tags.len() - 1) as usize {
                writeln!(output_stream, ",")?;
            } else {
                write!(output_stream, "\n")?;
            }
        }
        writeln!(output_stream, "\t\t]")?;
        writeln!(output_stream, "\t}},")?;

        writeln!(output_stream, "\t\"aln_tag_desc\" : {{")?;
        writeln!(
            output_stream,
            "\t\t\"label\" : \"{:?}\",",
            prelude.aln_tags.label
        )?;

        writeln!(output_stream, "\t\t\"tag_desc\" : [")?;
        for (i, td) in prelude.aln_tags.tags.iter().enumerate() {
            writeln!(
                output_stream,
                "\t\t\t{{\n\t\t\t \"name\" : \"{}\",",
                td.name
            )?;
            write!(
                output_stream,
                "\t\t\t \"desc\" : \"{:?}\"\n\t\t\t}}",
                td.typeid
            )?;
            if i < (prelude.aln_tags.tags.len() - 1) as usize {
                writeln!(output_stream, ",")?;
            } else {
                write!(output_stream, "\n")?;
            }
        }
        writeln!(output_stream, "\t\t]")?;
        writeln!(output_stream, "\t}}")?;

        //writeln!(output_stream, "{:?}", prelude.file_tags)?;
        /*
        pub file_tags: TagSection,
        pub read_tags: TagSection,
        pub aln_tags: TagSection,
        */
        writeln!(output_stream, "}},")?;
    }

    writeln!(output_stream, "\"mapped_records\" : [")?;
    let tag_context = prelude.get_record_context::<PiscemBulkRecordContext>()?;
    let num_chunks = 10; //prelude.hdr.num_chunks
    for chunk_num in 0..(num_chunks) {
        let chunk =
            libradicl::chunk::Chunk::<PiscemBulkReadRecord>::from_bytes(&mut ifile, &tag_context);
        let nreads = chunk.reads.len();
        for (rnum, r) in chunk.reads.iter().enumerate() {
            writeln!(output_stream, "{{")?;
            writeln!(
                output_stream,
                "\t\"frag_type\" : \"{:?}\",",
                libradicl::rad_types::MappingType::from_u8(r.frag_type)
            )?;
            writeln!(output_stream, "\t\"alns\" : [")?;

            for i in 0..(r.refs.len()) {
                write!(
                    output_stream,
                    "\t\t{{\"ref\": {}, \"dir\": \"{:?}\", \"pos\": {}, \"flen\": {} }}",
                    r.refs[i], r.dirs[i], r.positions[i], r.frag_lengths[i]
                )?;
                if i < r.refs.len() - 1 {
                    write!(output_stream, ",\n")?;
                } else {
                    write!(output_stream, "\n")?;
                }
            }

            writeln!(output_stream, "\t]")?;
            write!(output_stream, "}}")?;
            if (chunk_num == num_chunks - 1) && (rnum == nreads - 1) {
                write!(output_stream, "\n")?;
            } else {
                write!(output_stream, ",\n")?;
            }
        }
    }
    writeln!(output_stream, "]")?;

    writeln!(output_stream, "}}")?;

    Ok(())
}
