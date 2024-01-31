import sys
from lib import GffFile, longestPeps
from Bio import SeqIO


def longest_cds(gff_file, fasta_file, protein_seqs, cds_file, output_file):
    _, longest_trans_name = longestPeps.longestPeps(gff_file, fasta_file, protein_seqs, output_file)
    _, _, _, fake_transcript_gene_map = GffFile.readGff(gff_file)
    seqs = []
    for seq_record in SeqIO.parse(cds_file, "fasta"):
        if seq_record.id in longest_trans_name:
            seq_record.id = fake_transcript_gene_map[seq_record.id]
            seqs.append(seq_record)
        else:
            continue
    SeqIO.write(seqs, output_file, "fasta")


if __name__ == '__main__':
    longest_cds(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
