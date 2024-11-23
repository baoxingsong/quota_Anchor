from . import GffFile, FastaFile
import sys
import logging
import numpy as np
from argparse import ArgumentParser


# baoxing.song@pku-iaas.edu.cn
logger = logging.getLogger('main.longestPeps')

def longestPeps(gffFile, fastaFile, proteinSeqs, outputFile):
    logger.info(f"generate {outputFile} start.")
    longest_trans_name = []
    chromosome_gene_dict, chromosome_gene_list, geneName_toChr_dict, _ = GffFile.readGff(gffFile)
    chromosome_names, fastas = FastaFile.readFastaFile(fastaFile)
    pep_names, pep_fastas = FastaFile.readFastaFile(proteinSeqs)
    GffFile.update_sequence_information(fastas, chromosome_gene_dict)
    output = open(outputFile, 'w')
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
                if chromosome_gene_dict[chromosome_name][gene_name].transcripts[0].name in pep_fastas:
                    output.write(">" + gene_name + "\n")
                    output.write(pep_fastas[chromosome_gene_dict[chromosome_name][gene_name].transcripts[0].name].seq)
                    output.write("\n")
                    longest_trans_name.append(chromosome_gene_dict[chromosome_name][gene_name].transcripts[0].name)
    output.close()
    logger.info(f"generate {outputFile} done!")
    return chromosome_gene_dict, longest_trans_name


if __name__ == '__main__':
    parser = ArgumentParser(description='Prepare file for the strand and WGD aware syntenic gene identification function implemented in AnchorWave')
    parser.add_argument("-g", "--GFF",
                        dest="GffFile",
                        type=str,
                        default="",
                        help="Genome annotation in GFF formation")
    parser.add_argument("-f", "--fastaFile",
                        dest="fastaFile",
                        type=str,
                        default="",
                        help="Genome sequences in FASTA format")
    parser.add_argument("-p", "--proteinSeqs",
                        dest="proteinSeqs",
                        type=str,
                        default="",
                        help="Protein sequences in FASTA format")

    parser.add_argument("-o", "--output",
                        dest="outputFile",
                        type=str,
                        default="",
                        help="output file")

    args = parser.parse_args()

    if args.GffFile == "":
        print("Error: please specify --GffFile", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.fastaFile == "":
        print("Error: please specify --fastaFile", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.proteinSeqs == "":
        print("Error: please specify --proteinSeqs", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.outputFile == "":
        print("Error: please specify --output", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    gffFile = args.GffFile
    fastaFile = args.fastaFile
    proteinSeqs = args.proteinSeqs
    outputFile = args.outputFile

    longestPeps(gffFile, fastaFile, proteinSeqs, outputFile)
