# map_sRNA_reads
# Mapping sRNA reads from .fq to .fa
from Bio import SeqIO
from Bio.Seq import Seq

def read_fasta(file):
    reads = {}
    total_reads = 0
    for record in SeqIO.parse(file, "fasta"):
        header = record.id
        count = int(header.split('-')[-1])
        seq = str(record.seq)
        reads[seq] = count
        total_reads += count
    return reads, total_reads

def calculate_rpm(count, total_reads):
    return (float(count) / total_reads) * 1000000 if total_reads else 0

def map_sRNA_reads(transcript_seq, transcript_name, reads_dict, total_reads):
    results = []
    rev_seq = str(Seq(transcript_seq).reverse_complement())
    for read_seq, count in reads_dict.items():
        size = len(read_seq)

        pos = transcript_seq.find(read_seq)
        while pos != -1:
            start = pos + 1
            end = start + size - 1
            rpm = calculate_rpm(count, total_reads)
            results.append((read_seq, size, count, rpm, transcript_name, start, end, '+'))
            pos = transcript_seq.find(read_seq, pos + 1)

        pos_rev = rev_seq.find(read_seq)
        while pos_rev != -1:
            start = len(transcript_seq) - pos_rev - size + 1
            end = start + size - 1
            rpm = calculate_rpm(count, total_reads)
            results.append((read_seq, size, count, rpm, transcript_name, start, end, '-'))
            pos_rev = rev_seq.find(read_seq, pos_rev + 1)

    return results

if __name__ == "__main__":
    transcript_file = "NbTAS3-1.fa"
    collapsed_reads_file = "AC-HS-92_1_collapsed.fa"
    output_file = "mapping_summary_NbTAS3-1.tsv"

    transcript_record = SeqIO.read(transcript_file, "fasta")
    transcript_seq = str(transcript_record.seq)
    transcript_name = transcript_record.id

    reads_dict, total_reads = read_fasta(collapsed_reads_file)
    total_unique_reads = len(reads_dict)

    mappings = map_sRNA_reads(transcript_seq, transcript_name, reads_dict, total_reads)

    with open(output_file, "w") as out:
        out.write("{} total reads\n".format(total_reads))
        out.write("{} total different reads\n".format(total_unique_reads))
        out.write("smallRNA\tsize\treads\tRPMs\tgene\tstart\tend\tstrand\n")
        for entry in mappings:
            smallRNA, size, count, rpm, gene, start, end, strand = entry
            out.write("{}\t{}\t{}\t{:.6f}\t{}\t{}\t{}\t{}\n".format(
                smallRNA, size, count, rpm, gene, start, end, strand))

    print("Mapping complete. Results saved to {}".format(output_file))
