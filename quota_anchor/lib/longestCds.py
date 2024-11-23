import sys
import logging
from . import GffFile, FastaFile
from Bio import SeqIO
import numpy as np


logger = logging.getLogger('main.longestCds')
def longestCds(gffFile, fastaFile):
    longest_trans_name = []
    chromosome_gene_dict, _, _, _ = GffFile.readGff(gffFile)
    chromosome_names, fastas = FastaFile.readFastaFile(fastaFile)
    GffFile.update_sequence_information(fastas, chromosome_gene_dict)

    # get the gene with longest length
    # delete ORF-shift transcript
    # delete gene without transcript
    # keep only one transcript for each gene
    for chromosome_name in chromosome_names:
        gene_names_to_delete = []
        if (chromosome_name in chromosome_gene_dict) and len(chromosome_gene_dict[chromosome_name])>1:
            for gene_name in chromosome_gene_dict[chromosome_name]:
                if len(chromosome_gene_dict[chromosome_name][gene_name].transcripts) > 1:
                    longest = 0
                    for transcript in chromosome_gene_dict[chromosome_name][gene_name].transcripts:
                        if longest < len(transcript.cds_sequence):
                            longest = len(transcript.cds_sequence)

                    transcript_number = 0
                    for transcript in chromosome_gene_dict[chromosome_name][gene_name].transcripts:
                        if longest > len(transcript.cds_sequence):
                            chromosome_gene_dict[chromosome_name][gene_name].transcripts = np.delete(chromosome_gene_dict[chromosome_name][gene_name].transcripts, transcript_number, 0)
                            transcript_number = transcript_number - 1
                        transcript_number = transcript_number + 1

                # delete transcript so that the there is only one transcript left for each gene
                while len(chromosome_gene_dict[chromosome_name][gene_name].transcripts) > 1:
                    chromosome_gene_dict[chromosome_name][gene_name].transcripts = np.delete(chromosome_gene_dict[chromosome_name][gene_name].transcripts, 1, 0)
                # delete those genes contrains no transcript
                if len(chromosome_gene_dict[chromosome_name][gene_name].transcripts) == 0:
                    gene_names_to_delete.append(gene_name)

            for gene_name in gene_names_to_delete:
                del chromosome_gene_dict[chromosome_name][gene_name]

    for chromosome_name in chromosome_names:
        if (chromosome_name in chromosome_gene_dict) and len(chromosome_gene_dict[chromosome_name])>1:
            for gene_name in chromosome_gene_dict[chromosome_name]:
                longest_trans_name.append(chromosome_gene_dict[chromosome_name][gene_name].transcripts[0].name)

    return chromosome_gene_dict, longest_trans_name

# TODO: accelerate
def longest_cds(gff_file, fasta_file, cds_file, output_file):
    logger.info(f"generate {output_file} start.")
    output_file_handle = open(output_file, "a+")
    _, longest_trans_name = longestCds(gff_file, fasta_file)
    longest_trans_name = set(longest_trans_name)
    _, _, _, fake_transcript_gene_map = GffFile.readGff(gff_file)
    for seq_record in SeqIO.parse(cds_file, "fasta"):
        if seq_record.id in longest_trans_name:
            seq_record.id = fake_transcript_gene_map[seq_record.id]
            SeqIO.write(seq_record, output_file_handle, "fasta")
        else:
            continue
    logger.info(f"generate {output_file} done!")

if __name__ == '__main__':
    longest_cds(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
