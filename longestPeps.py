#!python

import NucleotideCodeSubstitution
import FastaFile
import GffFile
import MyUtil
import sys
import numpy as np

import math

# song@mpipz.mpg.de

_buckets = []

def read_data(gffFile, fastaFile, proteinSeqs, outputFile):
    chromosome_gene_dict, chromosome_gene_list, geneName_toChr_dict = GffFile.readGff( gffFile )
    chromosome_names, fastas = FastaFile.readFastaFile( fastaFile )
    pep_names, pep_fastas = FastaFile.readFastaFile(proteinSeqs)
    GffFile.update_sequence_information(fastas, chromosome_gene_dict)
    output = open(outputFile, 'w')
    #get the gene with longest length
    #delete ORF-shift transcript
    #delete gene without transcript
    #keep only one transcript for each gene
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

                ##delete transcript so that the there is only one transcript left for each gene
                while len(chromosome_gene_dict[chromosome_name][gene_name].transcripts)>1:
                    chromosome_gene_dict[chromosome_name][gene_name].transcripts = np.delete(chromosome_gene_dict[chromosome_name][gene_name].transcripts, 1, 0)
                #delete those genes contrains no transcript
                if len(chromosome_gene_dict[chromosome_name][gene_name].transcripts) == 0 :
                    gene_names_to_delete.append(gene_name)

            for gene_name in gene_names_to_delete:
                del chromosome_gene_dict[chromosome_name][gene_name]

    for chromosome_name in chromosome_names:
        if (chromosome_name in chromosome_gene_dict) and len(chromosome_gene_dict[chromosome_name])>1:
            for gene_name in chromosome_gene_dict[chromosome_name]:
                output.write(">" + gene_name + "\n")
                output.write(pep_fastas[chromosome_gene_dict[chromosome_name][gene_name].transcripts[0].name].seq)
                output.write("\n")

    output.close()
    return chromosome_gene_dict

print ("begin to run")
read_data(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
print ("stop running")
