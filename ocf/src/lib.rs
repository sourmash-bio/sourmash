/*
Copyright (c) 2018 Pierre Marijon <pierre.marijon@inria.fr>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Originally from https://github.com/natir/yacrd/blob/3fc6ef8b5b51256f0c4bc45b8056167acf34fa58/src/file.rs
Changes:
  - make bzip2 and lzma support optional
*/

/* crates use */
use cfg_if::cfg_if;
use enum_primitive::{
    enum_from_primitive, enum_from_primitive_impl, enum_from_primitive_impl_ty, FromPrimitive,
};
use failure::{Error, Fail};
use flate2;

/* standard use */
use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter};

enum_from_primitive! {
    #[repr(u64)]
    #[derive(Debug, PartialEq)]
    pub enum CompressionFormat {
        Gzip = 0x1F8B,
        Bzip = 0x425A,
        Lzma = 0xFD377A585A,
        No,
    }
}

#[derive(Debug, Fail)]
pub enum OCFError {
    #[fail(display = "Feature disabled, enabled it during compilation")]
    FeatureDisabled,
}

pub fn get_input(input_name: &str) -> Result<(Box<dyn io::Read>, CompressionFormat), Error> {
    // choose std::io::stdin or open file
    if input_name == "-" {
        Ok((Box::new(get_readable(input_name)), CompressionFormat::No))
    } else {
        get_readable_file(input_name)
    }
}

pub fn get_readable_file(
    input_name: &str,
) -> Result<(Box<dyn io::Read>, CompressionFormat), Error> {
    let raw_input = get_readable(input_name);

    // check compression
    let compression = get_compression(raw_input);

    // return readable and compression status
    match compression {
        CompressionFormat::Gzip => Ok((
            Box::new(flate2::read::GzDecoder::new(get_readable(input_name))),
            CompressionFormat::Gzip,
        )),
        CompressionFormat::Bzip => new_bz2_decoder(get_readable(input_name)),
        CompressionFormat::Lzma => new_lzma_decoder(get_readable(input_name)),
        CompressionFormat::No => Ok((Box::new(get_readable(input_name)), CompressionFormat::No)),
    }
}

pub fn get_readable(input_name: &str) -> Box<dyn io::Read> {
    match input_name {
        "-" => Box::new(BufReader::new(io::stdin())),
        _ => Box::new(BufReader::new(
            File::open(input_name).expect(&format!("Can't open input file {}", input_name)),
        )),
    }
}

fn get_compression(mut in_stream: Box<dyn io::Read>) -> CompressionFormat {
    let mut buf = vec![0u8; 5];

    in_stream
        .read_exact(&mut buf)
        .expect("Error durring reading first bit of file");

    let mut five_bit_val: u64 = 0;
    for i in 0..5 {
        five_bit_val |= (buf[i] as u64) << 8 * (4 - i);
    }
    if CompressionFormat::from_u64(five_bit_val) == Some(CompressionFormat::Lzma) {
        return CompressionFormat::Lzma;
    }

    let mut two_bit_val: u64 = 0;
    for i in 0..2 {
        two_bit_val |= (buf[i] as u64) << 8 * (1 - i);
    }

    match CompressionFormat::from_u64(two_bit_val) {
        e @ Some(CompressionFormat::Gzip) | e @ Some(CompressionFormat::Bzip) => e.unwrap(),
        _ => CompressionFormat::No,
    }
}

cfg_if! {
    if #[cfg(feature = "bz2")] {
        use bzip2;

        fn new_bz2_encoder(out: Box<dyn io::Write>) -> Result<Box<dyn io::Write>, Error> {
            Ok(Box::new(bzip2::write::BzEncoder::new(
                out,
                bzip2::Compression::Best,
            )))
        }

        fn new_bz2_decoder(
            inp: Box<dyn io::Read>,
        ) -> Result<(Box<dyn io::Read>, CompressionFormat), Error> {
            use bzip2;
            Ok((
                Box::new(bzip2::read::BzDecoder::new(inp)),
                CompressionFormat::Bzip,
            ))
        }
    } else {
        fn new_bz2_encoder(_: Box<dyn io::Write>) -> Result<Box<dyn io::Write>, Error> {
            Err(OCFError::FeatureDisabled.into())
        }

        fn new_bz2_decoder(_: Box<dyn io::Read>) -> Result<(Box<dyn io::Read>, CompressionFormat), Error> {
            Err(OCFError::FeatureDisabled.into())
        }
    }
}

cfg_if! {
    if #[cfg(feature = "lzma")] {
      use xz2;

      fn new_lzma_encoder(out: Box<dyn io::Write>) -> Result<Box<dyn io::Write>, Error> {
          Ok(Box::new(xz2::write::XzEncoder::new(out, 9)))
      }

      fn new_lzma_decoder(
          inp: Box<dyn io::Read>,
      ) -> Result<(Box<dyn io::Read>, CompressionFormat), Error> {
          use xz2;
          Ok((
              Box::new(xz2::read::XzDecoder::new(inp)),
              CompressionFormat::Lzma,
          ))
      }
    } else {
      fn new_lzma_encoder(_: Box<dyn io::Write>) -> Result<Box<dyn io::Write>, Error> {
          Err(OCFError::FeatureDisabled.into())
      }

      fn new_lzma_decoder(_: Box<dyn io::Read>) -> Result<(Box<dyn io::Read>, CompressionFormat), Error> {
          Err(OCFError::FeatureDisabled.into())
      }
    }
}

pub fn get_output(
    output_name: &str,
    format: CompressionFormat,
) -> Result<Box<dyn io::Write>, Error> {
    match format {
        CompressionFormat::Gzip => Ok(Box::new(flate2::write::GzEncoder::new(
            get_writable(output_name),
            flate2::Compression::best(),
        ))),
        CompressionFormat::Bzip => new_bz2_encoder(get_writable(output_name)),
        CompressionFormat::Lzma => new_lzma_encoder(get_writable(output_name)),
        CompressionFormat::No => Ok(Box::new(get_writable(output_name))),
    }
}

pub fn choose_compression(
    input_compression: CompressionFormat,
    compression_set: bool,
    compression_value: &str,
) -> CompressionFormat {
    if !compression_set {
        return input_compression;
    }

    match compression_value {
        "gzip" => CompressionFormat::Gzip,
        "bzip2" => CompressionFormat::Bzip,
        "lzma" => CompressionFormat::Lzma,
        _ => CompressionFormat::No,
    }
}

fn get_writable(output_name: &str) -> Box<dyn io::Write> {
    match output_name {
        "-" => Box::new(BufWriter::new(io::stdout())),
        _ => Box::new(BufWriter::new(
            File::create(output_name).expect(&format!("Can't open output file {}", output_name)),
        )),
    }
}

#[cfg(test)]
mod test {

    use super::*;

    const GZIP_FILE: &'static [u8] = &[0o037, 0o213, 0o0, 0o0, 0o0];
    const BZIP_FILE: &'static [u8] = &[0o102, 0o132, 0o0, 0o0, 0o0];
    const LZMA_FILE: &'static [u8] = &[0o375, 0o067, 0o172, 0o130, 0o132];

    #[test]
    fn compression_from_file() {
        assert_eq!(
            get_compression(Box::new(GZIP_FILE)),
            CompressionFormat::Gzip
        );
        assert_eq!(
            get_compression(Box::new(BZIP_FILE)),
            CompressionFormat::Bzip
        );
        assert_eq!(
            get_compression(Box::new(LZMA_FILE)),
            CompressionFormat::Lzma
        );
    }

    #[test]
    fn compression_from_input_or_cli() {
        assert_eq!(
            choose_compression(CompressionFormat::Gzip, false, "_"),
            CompressionFormat::Gzip
        );
        assert_eq!(
            choose_compression(CompressionFormat::Bzip, false, "_"),
            CompressionFormat::Bzip
        );
        assert_eq!(
            choose_compression(CompressionFormat::Lzma, false, "_"),
            CompressionFormat::Lzma
        );
        assert_eq!(
            choose_compression(CompressionFormat::No, true, "gzip"),
            CompressionFormat::Gzip
        );
        assert_eq!(
            choose_compression(CompressionFormat::No, true, "bzip2"),
            CompressionFormat::Bzip
        );
        assert_eq!(
            choose_compression(CompressionFormat::No, true, "lzma"),
            CompressionFormat::Lzma
        );
    }
}
