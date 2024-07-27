use biocirc::sequence::BioType;
use biocirc::FastaReader;

#[test]
fn read_fasta() {
    let fastareader = FastaReader {
        file: "D:\\BioCirc\\tests\\data\\example.fasta".to_string(),
        biotype: BioType::Dna,
    };
    let sequences = fastareader.read().unwrap();
    
    assert_eq!(sequences.len(), 3);
    println!("{:?}", sequences[0].sequence);
    println!("{:?}", sequences[0].id);
    println!("{:?}", sequences[0].description);
    println!("{:?}", sequences[1].sequence);
    println!("{:?}", sequences[1].id);
    println!("{:?}", sequences[1].description);
}
