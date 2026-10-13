//! NEXUS lexical helpers shared by the report writer and reader: token quoting, a
//! quote-aware line tokenizer, and the plain / zstd file sinks and sources.

use anyhow::{bail, Context, Result};
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/// First bytes of a zstd frame.
pub const ZSTD_MAGIC: [u8; 4] = [0x28, 0xB5, 0x2F, 0xFD];

/// Characters a token may hold unquoted. Everything else (whitespace, NEXUS
/// punctuation, `_`, which NEXUS readers turn into a space, `|`, `:` ...) is quoted.
fn is_bare(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '.' || c == '-' || c == '+'
}

/// Quote `s` for a NEXUS file when needed: `'...'` with every `'` doubled. Control
/// characters (which cannot occur in a FASTA header or FASTQ id, but would break the
/// line-oriented reader) become spaces.
pub(crate) fn quote(s: &str) -> Cow<'_, str> {
    if !s.is_empty() && s.chars().all(is_bare) {
        return Cow::Borrowed(s);
    }
    let mut out = String::with_capacity(s.len() + 2);
    out.push('\'');
    for c in s.chars() {
        match c {
            '\'' => out.push_str("''"),
            c if c.is_control() => out.push(' '),
            c => out.push(c),
        }
    }
    out.push('\'');
    Cow::Owned(out)
}

/// Split one line into tokens: whitespace-separated words, `'quoted'` strings (with
/// `''` for a quote), and `[comments]` outside quotes dropped. A trailing `;` on a bare
/// word is split off as its own token.
pub(crate) fn tokens(line: &str) -> Result<Vec<String>> {
    let mut out = Vec::new();
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        if c.is_whitespace() {
            continue;
        }
        match c {
            '\'' => {
                let mut tok = String::new();
                loop {
                    match chars.next() {
                        Some('\'') => {
                            if chars.peek() == Some(&'\'') {
                                chars.next();
                                tok.push('\'');
                            } else {
                                break;
                            }
                        }
                        Some(ch) => tok.push(ch),
                        None => bail!("unterminated quote in: {}", line),
                    }
                }
                out.push(tok);
            }
            '[' => {
                let mut depth = 1;
                for ch in chars.by_ref() {
                    match ch {
                        '[' => depth += 1,
                        ']' => {
                            depth -= 1;
                            if depth == 0 {
                                break;
                            }
                        }
                        _ => {}
                    }
                }
                if depth != 0 {
                    bail!("unterminated comment in: {}", line);
                }
            }
            ';' | ',' | '=' => out.push(c.to_string()),
            _ => {
                let mut tok = String::new();
                tok.push(c);
                while let Some(&n) = chars.peek() {
                    if n.is_whitespace()
                        || n == ';'
                        || n == '\''
                        || n == '['
                        || n == ','
                        || n == '='
                    {
                        break;
                    }
                    tok.push(n);
                    chars.next();
                }
                out.push(tok);
            }
        }
    }
    Ok(out)
}

/// Open a report for reading, transparently decoding zstd when the file starts with
/// its magic bytes.
pub fn open_report(path: &str) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("failed to open report '{}'", path))?;
    let mut reader = BufReader::new(file);
    let head = reader.fill_buf()?;
    if head.starts_with(&ZSTD_MAGIC) {
        let dec = zstd::Decoder::with_buffer(reader)
            .with_context(|| format!("failed to start zstd decoding of '{}'", path))?;
        Ok(Box::new(BufReader::new(dec)))
    } else {
        Ok(Box::new(reader))
    }
}

/// An output file, plain or zstd-compressed. Must be [`Sink::finish`]ed: dropping a
/// zstd encoder loses the final frame.
#[must_use = "a Sink must be finished, or the file stays incomplete"]
pub enum Sink {
    Plain(BufWriter<File>),
    Zstd(zstd::Encoder<'static, BufWriter<File>>),
}

impl Sink {
    pub fn create(path: &str, compress: bool) -> Result<Self> {
        let file = File::create(path).with_context(|| format!("failed to create '{}'", path))?;
        let buf = BufWriter::new(file);
        if compress {
            let enc = zstd::Encoder::new(buf, 0)
                .with_context(|| format!("failed to start zstd encoding of '{}'", path))?;
            Ok(Self::Zstd(enc))
        } else {
            Ok(Self::Plain(buf))
        }
    }

    /// Flush everything to disk.
    pub fn finish(self) -> Result<()> {
        match self {
            Self::Plain(mut w) => w.flush()?,
            Self::Zstd(enc) => enc.finish()?.flush()?,
        }
        Ok(())
    }
}

impl Write for Sink {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(w) => w.write(buf),
            Self::Zstd(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Plain(w) => w.flush(),
            Self::Zstd(w) => w.flush(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn quote_round_trips_awkward_labels() {
        for s in [
            "a b",
            "x|y:z",
            "it's",
            "_u_",
            "[1]",
            "",
            "a;b",
            "NC_000913.3",
            "plain",
            "1e-5",
        ] {
            let q = quote(s);
            let toks = tokens(&q).unwrap();
            assert_eq!(toks, vec![s.to_string()], "quoting {s:?} gave {q}");
        }
        assert_eq!(quote("plain"), "plain");
        assert_eq!(quote("a b"), "'a b'");
        assert_eq!(quote("it's"), "'it''s'");
        assert_eq!(quote(""), "''");
        assert_eq!(quote("a\nb"), "'a b'");
    }

    #[test]
    fn tokens_split_matrix_rows_and_drop_comments() {
        let t = tokens("  7 'NC 1' 3.5e-1 [note] 12 ;").unwrap();
        assert_eq!(t, vec!["7", "NC 1", "3.5e-1", "12", ";"]);
        let t = tokens("DIMENSIONS NTAX=3;").unwrap();
        assert_eq!(t, vec!["DIMENSIONS", "NTAX", "=", "3", ";"]);
        let t = tokens("TRANSLATE 0 root, 1 '5 62';").unwrap();
        assert_eq!(t, vec!["TRANSLATE", "0", "root", ",", "1", "5 62", ";"]);
        assert!(tokens("'open").is_err());
        assert!(tokens("[open").is_err());
    }

    #[test]
    fn sink_round_trips_plain_and_zstd() {
        let dir = tempfile::tempdir().unwrap();
        for compress in [false, true] {
            let path = dir
                .path()
                .join(if compress { "r.zst" } else { "r.txt" })
                .to_str()
                .unwrap()
                .to_string();
            let mut sink = Sink::create(&path, compress).unwrap();
            writeln!(sink, "#NEXUS").unwrap();
            writeln!(sink, "line two").unwrap();
            sink.finish().unwrap();
            let raw = std::fs::read(&path).unwrap();
            assert_eq!(raw.starts_with(&ZSTD_MAGIC), compress);
            let mut text = String::new();
            open_report(&path)
                .unwrap()
                .read_to_string(&mut text)
                .unwrap();
            assert_eq!(text, "#NEXUS\nline two\n");
        }
    }
}
