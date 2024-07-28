pub mod codon;
pub mod sequence;

use sequence::{BioType, Sequence};
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

/// 接下来是对于生物序列对象的封装
/// SeqRecode
/// biopython 应该没有版权法吧
/// 包括： 序列名称、类型、描述、序列

/// 用于封装 fasta 对象数据
/// 结构包括：id、name、description、sequence
pub struct FastaRecode {
    pub id: String,
    pub description: String,
    pub sequence: Sequence,
}

impl FastaRecode {
    /// 新建 FastaRecode 对象
    pub fn new(id: String, description: String, sequence: Sequence) -> Self {
        FastaRecode {
            id,
            description,
            sequence,
        }
    }
}

pub struct FastaReader {
    pub file: String,
    pub biotype: BioType,
}

impl FastaReader {
    pub fn read(&self) -> io::Result<Vec<FastaRecode>> {
        let path = self.file.clone();
        let file = File::open(&path)?;
        let reader = io::BufReader::new(file);
        let mut sequence_vec: Vec<FastaRecode> = Vec::new();
        let mut fasta_recode: Option<FastaRecode> = None;

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('>') {
                if let Some(rec) = fasta_recode.take() {
                    sequence_vec.push(rec);
                }
                let (id, description) = line[1..].split_once(' ').unwrap_or((&line[1..], ""));
                fasta_recode = Some(FastaRecode::new(
                    id.to_string(),
                    description.to_string(),
                    Sequence::new(self.biotype.clone(), String::new()),
                ));
            } else {
                if let Some(ref mut rec) = fasta_recode {
                    rec.sequence.seq.push_str(&line);
                }
            }
        }

        // Push the last record if exists
        if let Some(rec) = fasta_recode {
            sequence_vec.push(rec);
        }

        Ok(sequence_vec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(test)]

    fn base_function() {
        let result: Sequence = Sequence::new(BioType::Dna, String::from("ACGTTCGA"));
        assert_eq!(result.len(), 8);
        assert_eq!(result.index(2), 'G');
    }

    #[cfg(test)]

    fn string() {
        let mut sequence1: Sequence = Sequence::new(BioType::Dna, String::from("AGCT"));
        println!("Sequence: {}", sequence1);
        println!("index 2: {}", sequence1.index(2));
        sequence1.change(2, 'R');
        println!("change after: {}", sequence1);

        let sequence2: Sequence = Sequence::new(BioType::Dna, String::from("TCGA"));
        let sequence3: Sequence = Sequence::new(BioType::Dna, String::from("AGRTTCGA"));
        assert_eq!(sequence1 + sequence2, sequence3);
        let mut sequence1: Sequence = Sequence::new(BioType::Dna, String::from("AGCT"));
        sequence1.push('G');
        assert_eq!(
            sequence1,
            Sequence::new(BioType::Dna, String::from("AGCTG"))
        );
    }

    #[cfg(test)]

    fn trans() {
        let sequence1: Sequence = Sequence::new(BioType::Dna, String::from("AGCTTCGAA"));

        println!("{:?}", sequence1.transcribe());
        println!("{:?}", sequence1.translate());
        println!("{:?}", sequence1.reverse_complementary());
        println!("{:?}", sequence1.complementary());
    }

    #[cfg(test)]

    fn seqrecode() {
        let mut recode1 = FastaRecode {
            id: String::from("recode1"),
            description: String::new(),
            sequence: Sequence::new(BioType::Rna, String::from("AGCUTUCG")),
        };

        let recorde2 = FastaRecode::new(
            String::from("recode2"),
            String::new(),
            Sequence::new(BioType::Rna, String::from("AGCUTUCGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")),
        );
        recode1.sequence = Sequence::new(BioType::Dna, String::from("AGCTCGT"));
        println!("{}", recorde2.sequence);
    }
}
