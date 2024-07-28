use criterion::{black_box, criterion_group, criterion_main, Criterion};
use biocirc::sequence::BioType;
use biocirc::FastaReader;


fn read_fasta() {
    let fastareader = FastaReader {
        file: "tests\\data\\example.fasta".to_string(),
        biotype: BioType::Dna,
    };
    let sequences = fastareader.read().unwrap();
    assert_eq!(sequences.len(), 3);
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("fasta open", |b| b.iter(|| read_fasta()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);