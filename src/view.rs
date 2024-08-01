use anyhow::bail;
use clap::{Parser, ValueEnum};
use libradicl::record::{
    AlevinFryReadRecord, AlevinFryRecordContext, PiscemBulkReadRecord, PiscemBulkRecordContext,
};
use needletail::bitkmer::*;
use std::io;
use std::io::{BufReader, Write};
use tracing::error;

/// The types of RAD files supported
#[derive(Clone, Debug, PartialEq, ValueEnum)]
pub enum RadFileType {
    Bulk,
    SingleCell,
    Unknown,
}

/// options related to printing a RAD file
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct ViewOpts {
    /// the input RAD file to print
    #[arg(short, long, required = true)]
    input: std::path::PathBuf,

    /// output file where the JSON format RAD file will be written;
    /// if not provided, the output will be written to standard out.
    #[arg(short, long)]
    output: Option<std::path::PathBuf>,

    /// the type of input RAD file
    #[arg(short, long)]
    rad_type: RadFileType,

    /// use the reference name rather than ID in the mapped records
    #[arg(long)]
    use_ref_name: bool,

    /// skip printing the header and file-level tags (i.e. only print the mapping reacords)
    #[arg(long)]
    no_header: bool,

    /// print the records from at most this many chunks
    #[arg(long)]
    max_chunks: Option<usize>,
}

/// **NOTE**: This representation is a hack and we should think of
/// a better way to handle generic information over these.
pub struct ExtraRecordInfo<'a> {
    pub bc_len: usize,
    pub umi_len: usize,
    pub use_ref_name: bool,
    pub prelude: &'a libradicl::header::RadPrelude,
    pub max_chunks: Option<usize>,
}

impl<'a> ExtraRecordInfo<'a> {
    /// Provides the ability to use the header to lookup
    /// the name of a target given its ID.
    pub fn ref_name(&self, i: usize) -> &str {
        &self.prelude.hdr.ref_names[i]
    }
}

/// The ability to write mapping records of different types
pub trait WriteMappingRecord {
    fn write_records(
        &self,
        ctx: &ExtraRecordInfo,
        output_stream: &mut Box<dyn Write>,
    ) -> anyhow::Result<()>;
}

impl WriteMappingRecord for libradicl::record::PiscemBulkReadRecord {
    fn write_records(
        &self,
        ctx: &ExtraRecordInfo,
        output_stream: &mut Box<dyn Write>,
    ) -> anyhow::Result<()> {
        writeln!(output_stream, "{{")?;
        writeln!(
            output_stream,
            " \"frag_type\" : \"{:?}\",",
            libradicl::rad_types::MappingType::from_u8(self.frag_type)
        )?;
        writeln!(output_stream, " \"alns\" : [")?;

        for i in 0..(self.refs.len()) {
            if ctx.use_ref_name {
                write!(
                    output_stream,
                    "  {{\"ref\": {:?}, \"dir\": \"{:?}\", \"pos\": {}, \"flen\": {} }}",
                    ctx.ref_name(self.refs[i] as usize),
                    self.dirs[i],
                    self.positions[i],
                    self.frag_lengths[i]
                )?;
            } else {
                write!(
                    output_stream,
                    "  {{\"ref\": {}, \"dir\": \"{:?}\", \"pos\": {}, \"flen\": {} }}",
                    self.refs[i], self.dirs[i], self.positions[i], self.frag_lengths[i]
                )?;
            }

            if i < self.refs.len() - 1 {
                writeln!(output_stream, ",")?;
            } else {
                writeln!(output_stream)?;
            }
        }

        writeln!(output_stream, " ]")?;
        write!(output_stream, "}}")?;
        Ok(())
    }
}

impl WriteMappingRecord for libradicl::record::AlevinFryReadRecord {
    fn write_records(
        &self,
        ctx: &ExtraRecordInfo,
        output_stream: &mut Box<dyn Write>,
    ) -> anyhow::Result<()> {
        let bc_mer: BitKmer = (self.bc, ctx.bc_len as u8);
        let umi_mer: BitKmer = (self.umi, ctx.umi_len as u8);

        writeln!(output_stream, "{{")?;
        writeln!(
            output_stream,
            " \"barcode\" : {:?},\n \"umi\" : {:?},",
            unsafe { std::str::from_utf8_unchecked(&bitmer_to_bytes(bc_mer)[..]) },
            unsafe { std::str::from_utf8_unchecked(&bitmer_to_bytes(umi_mer)[..]) },
        )?;
        writeln!(output_stream, " \"alns\" : [")?;

        for i in 0..(self.refs.len()) {
            if ctx.use_ref_name {
                write!(
                    output_stream,
                    "  {{\"ref\": {:?}, \"dir\": \"{}\" }}",
                    ctx.ref_name(self.refs[i] as usize),
                    if self.dirs[i] { "fw" } else { "rc" }
                )?;
            } else {
                write!(
                    output_stream,
                    "  {{\"ref\": {}, \"dir\": \"{}\" }}",
                    self.refs[i],
                    if self.dirs[i] { "fw" } else { "rc" }
                )?;
            }

            if i < self.refs.len() - 1 {
                writeln!(output_stream, ",")?;
            } else {
                writeln!(output_stream)?;
            }
        }

        writeln!(output_stream, " ]")?;
        write!(output_stream, "}}")?;
        Ok(())
    }
}

pub fn write_records<
    RecordContext: std::fmt::Debug + Clone + libradicl::record::RecordContext,
    RecordType: std::fmt::Debug
        + libradicl::record::MappedRecord<ParsingContext = RecordContext>
        + WriteMappingRecord,
    R: std::io::BufRead,
>(
    prelude: &libradicl::header::RadPrelude,
    extra_record_info: &ExtraRecordInfo,
    ifile: &mut R,
    output_stream: &mut Box<dyn Write>,
) -> anyhow::Result<()> {
    let tag_context = prelude.get_record_context::<RecordContext>()?;
    let total_chunks = if prelude.hdr.num_chunks > 0 {
        prelude.hdr.num_chunks as usize
    } else {
        usize::MAX - 1
    };
    let mut chunk_num = 0;

    // output either the requested number of chunks, or all of them if no
    // request is provided (but never try to print more than the total).
    let num_chunks = extra_record_info
        .max_chunks
        .unwrap_or(total_chunks)
        .min(total_chunks);

    while chunk_num < num_chunks
        && libradicl::utils::has_data_left(ifile).expect("encountered error reading input file")
    {
        // write out each chunk.
        let chunk = libradicl::chunk::Chunk::<RecordType>::from_bytes(ifile, &tag_context);
        let nreads = chunk.reads.len();
        for (rnum, r) in chunk.reads.iter().enumerate() {
            r.write_records(extra_record_info, output_stream)?;
            if (chunk_num == num_chunks - 1) && (rnum == nreads - 1) {
                writeln!(output_stream)?;
            } else {
                writeln!(output_stream, ",")?;
            }
        }
        chunk_num += 1;
    }
    Ok(())
}

pub fn write_header(
    prelude: &libradicl::header::RadPrelude,
    file_tag_map: &libradicl::rad_types::TagMap,
    output_stream: &mut Box<dyn Write>,
) -> anyhow::Result<()> {
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

    writeln!(output_stream, " \"file_tag_desc\" : {{")?;
    writeln!(
        output_stream,
        "  \"label\" : \"{:?}\",",
        prelude.file_tags.label
    )?;

    writeln!(output_stream, "  \"tag_desc\" : [")?;
    for (i, td) in prelude.file_tags.tags.iter().enumerate() {
        writeln!(output_stream, "   {{\n    \"name\" : \"{}\",", td.name)?;
        write!(output_stream, "    \"desc\" : \"{:?}\"\n   }}", td.typeid)?;
        if i < (prelude.file_tags.tags.len() - 1) {
            writeln!(output_stream, ",")?;
        } else {
            writeln!(output_stream)?;
        }
    }
    writeln!(output_stream, "  ]")?;
    writeln!(output_stream, " }},")?;

    writeln!(output_stream, " \"read_tag_desc\" : {{")?;
    writeln!(
        output_stream,
        "  \"label\" : \"{:?}\",",
        prelude.read_tags.label
    )?;

    writeln!(output_stream, "  \"tag_desc\" : [")?;
    for (i, td) in prelude.read_tags.tags.iter().enumerate() {
        writeln!(output_stream, "   {{\n    \"name\" : \"{}\",", td.name)?;
        write!(output_stream, "    \"desc\" : \"{:?}\"\n   }}", td.typeid)?;
        if i < (prelude.read_tags.tags.len() - 1) {
            writeln!(output_stream, ",")?;
        } else {
            writeln!(output_stream)?;
        }
    }
    writeln!(output_stream, "  ]")?;
    writeln!(output_stream, " }},")?;

    writeln!(output_stream, " \"aln_tag_desc\" : {{")?;
    writeln!(
        output_stream,
        "  \"label\" : \"{:?}\",",
        prelude.aln_tags.label
    )?;

    writeln!(output_stream, "  \"tag_desc\" : [")?;
    for (i, td) in prelude.aln_tags.tags.iter().enumerate() {
        writeln!(output_stream, "   {{\n    \"name\" : \"{}\",", td.name)?;
        write!(output_stream, "    \"desc\" : \"{:?}\"\n   }}", td.typeid)?;
        if i < (prelude.aln_tags.tags.len() - 1) {
            writeln!(output_stream, ",")?;
        } else {
            writeln!(output_stream)?;
        }
    }
    writeln!(output_stream, "  ]")?;
    writeln!(output_stream, " }}")?;

    //writeln!(output_stream, "{:?}", prelude.file_tags)?;
    /*
    pub file_tags: TagSection,
    pub read_tags: TagSection,
    pub aln_tags: TagSection,
    */
    writeln!(output_stream, "}},")?;

    // file tags
    writeln!(output_stream, "\"file_tags\" : [")?;

    let nft = prelude.file_tags.tags.len();
    for (i, td) in prelude.file_tags.tags.iter().enumerate() {
        if let Some(tv) = file_tag_map.get(&td.name) {
            write!(
                output_stream,
                "{{ \"name\" : \"{}\", \"val\" : \"{:?}\" }}",
                &td.name, tv
            )?;
            if i == nft - 1 {
                writeln!(output_stream)?;
            } else {
                writeln!(output_stream, ",")?;
            }
        }
    }
    writeln!(output_stream, "],")?;
    Ok(())
}

pub fn view(view_opts: &ViewOpts) -> anyhow::Result<()> {
    if view_opts.rad_type == RadFileType::Unknown {
        error!("Unknown file type not yet supported");
        bail!("Unknown file type not yet supported");
    }

    let mut output_stream: Box<dyn Write> = match view_opts.output {
        Some(ref path) => std::fs::OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(path)
            .map(|f| Box::new(std::io::BufWriter::new(f)) as Box<dyn Write>)?,
        None => Box::new(io::stdout()),
    };

    let f = std::fs::File::open(&view_opts.input)?;
    let mut ifile = BufReader::new(f);
    let prelude = libradicl::header::RadPrelude::from_bytes(&mut ifile)?;
    let file_tag_map = prelude.file_tags.try_parse_tags_from_bytes(&mut ifile)?;

    writeln!(output_stream, "{{")?;

    if !view_opts.no_header {
        write_header(&prelude, &file_tag_map, &mut output_stream)?;
    }

    let mut extra_record_info = ExtraRecordInfo {
        bc_len: 0,
        umi_len: 0,
        use_ref_name: view_opts.use_ref_name,
        prelude: &prelude,
        max_chunks: view_opts.max_chunks,
    };

    writeln!(output_stream, "\"mapped_records\" : [")?;
    match view_opts.rad_type {
        RadFileType::Bulk => {
            write_records::<PiscemBulkRecordContext, PiscemBulkReadRecord, BufReader<std::fs::File>>(
                &prelude,
                &extra_record_info,
                &mut ifile,
                &mut output_stream,
            )?;
        }
        RadFileType::SingleCell => {
            let cblen: u64 = file_tag_map
                .get("cblen")
                .expect("tag map must contain \"cblen\" value")
                .try_into()?;

            let ulen: u64 = file_tag_map
                .get("ulen")
                .expect("tag map must contain \"ulen\" value")
                .try_into()?;

            extra_record_info.bc_len = cblen as usize;
            extra_record_info.umi_len = ulen as usize;

            write_records::<AlevinFryRecordContext, AlevinFryReadRecord, BufReader<std::fs::File>>(
                &prelude,
                &extra_record_info,
                &mut ifile,
                &mut output_stream,
            )?;
        }
        RadFileType::Unknown => bail!("Unknown RadFileType not supported yet"),
    }

    writeln!(output_stream, "]")?;
    writeln!(output_stream, "}}")?;

    Ok(())
}
