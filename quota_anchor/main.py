import configparser
import logging
import argparse
from pathlib import Path
import os
import sys
from .lib import base, collinearity, dotplot, circle, get_chr_length, line, classification_gene
from .lib import get_longest_pep, pre_collinearity, get_longest_cds, ks
from .kspeaks import kde, ks_fitting, trios, correct

base_dir = Path(__file__).resolve().parent


def run_get_longest_pep(parameter):
    global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)

    config_soft = configparser.ConfigParser()
    config_soft.read(os.path.join(base_dir, 'config_file/software_path.ini'))
    get_longest_pep.Longest(config_par,config_soft, parameter).run_all_process()

def run_get_longest_cds(parameter):
    global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)
    config_soft = configparser.ConfigParser()
    config_soft.read(os.path.join(base_dir, 'config_file/software_path.ini'))
    get_longest_cds.Longest(config_par,config_soft, parameter).run_all_process()

def run_pre_col(parameter):
    global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)
    config_soft = configparser.ConfigParser()
    config_soft.read(os.path.join(base_dir, 'config_file/software_path.ini'))
    pre_collinearity.Prepare(config_par,config_soft, parameter).run_all_process()


def run_col(parameter):
    global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)
     
    config_soft = configparser.ConfigParser()
    config_soft.read(os.path.join(base_dir, 'config_file/software_path.ini'))
    collinearity.Collinearity(config_par, config_soft, parameter).run()


def run_dotplot(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)
    dotplot.Dotplot(config_par, parameter).run()


def run_circle(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)
    circle.Circle(config_par, parameter).run()


def run_lens(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)
    get_chr_length.Lens(config_par, parameter).run()


def run_line(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)
    line.Line(config_par, parameter).run()


def run_ks(parameter):
    global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)
    config_soft = configparser.ConfigParser()
    config_soft.read(os.path.join(base_dir, 'config_file/software_path.ini'))
    ks.Ks(config_par,config_soft, parameter).first_layer_run()


def run_class_gene(parameter):
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)
    if parameter.unique:
        classification_gene.ClassGeneUnique(config_par, parameter).run()
    else:
        classification_gene.ClassGene(config_par, parameter).run()

def run_kde(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)
    kde.Kde(config_par, parameter).run()

def run_kf(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)
    ks_fitting.Kf(config_par, parameter).run()

def run_trios(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)
    trios.Trios(config_par, parameter).run()

def run_correct(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        base.file_empty(parameter.conf)
        config_par.read(parameter.conf)
    correct.Correct(config_par, parameter).run()

parser = argparse.ArgumentParser(description='Conduct strand and WGD aware syntenic gene identification for a pair of genomes using the longest path algorithm implemented in AnchorWave.', prog="quota_Anchor")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1rc')


subparsers = parser.add_subparsers(title='Gene collinearity analysis', dest='analysis')

# get the longest protein sequence file
parser_sub_longest_pep = subparsers.add_parser('longest_pep', help='Call gffread to generate the protein sequence of the species based on the genome and gff files. The longest transcripts are then extracted based on the gff file and the protein sequence file.', formatter_class=argparse.RawDescriptionHelpFormatter, description="""
    You can execute this command in three ways:
                                       
    1. Using configuration file:
       a)quota_Anchor longest_pep -c [\\?, example, help] >> longest_pep.conf 
       b)quota_Anchor longest_pep -c longest_pep.conf [--overwrite] [-merge merged.fa]
                                       
    2. Using command-line arguments:
       quota_Anchor longest_pep -f sorghum.fa,maize.fa -g sorghum.gff3,maize.gff3 
                                -p sb.p.fa,zm.p.fa -l sorghum.protein.fa,maize.protein.fa -t 2 
                                [--overwrite] [-merge merged.fa]
                                       
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
       
       quota_Anchor longest_pep -f sorghum.fa,maize.fa -g sorghum.gff3,maize.gff3 
                                -p sb.p.fa,zm.p.fa -l sorghum.protein.fa,maize.protein.fa -t 2 -c longest_pep.conf
                                [--overwrite] [-merge merged.fa]                                      
 """)
parser_sub_longest_pep.set_defaults(func=run_get_longest_pep)
parser_sub_longest_pep.add_argument('-c', '--conf', dest='conf', help="Configuration files have the lowest priority.", metavar="")
parser_sub_longest_pep.add_argument('-f', '--genome_file', dest='genome_file', help="Species genome file list(separator: ',').", metavar="")
parser_sub_longest_pep.add_argument('-g', '--gff_file', dest='gff_file', help="Species gff file list(separator: ',').", metavar="")
parser_sub_longest_pep.add_argument('-p', '--out_pep_file', dest='out_pep_file', help="Output species raw protein file list(separator: ',').", metavar="")
parser_sub_longest_pep.add_argument('-l', '--out_longest_pep_file', dest='out_longest_pep_file', help="Output species longest protein file list(separator: ',').", metavar="")
parser_sub_longest_pep.add_argument('-t', '--thread', dest='thread', help="Number of parallel processes.", type=int, metavar="")
parser_sub_longest_pep.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file.", action='store_true')
parser_sub_longest_pep.add_argument('-merge', '--merge', dest='merge', help="Optional, specify a file name, and merge all longest pep file content to this file.", metavar="")


# get the longest cds sequence file
parser_sub_longest_cds = subparsers.add_parser('longest_cds', help='Call gffread to generate the coding sequence of the species based on the genome and gff files. The longest cds are then extracted based on the gff file and the coding sequence file.', formatter_class=argparse.RawDescriptionHelpFormatter, description="""
    You can execute this command in three ways: 
                                       
    1. Using configuration file:
       a)quota_Anchor longest_cds -c [\\?, example, help] >> longest_cds.conf 
       b)quota_Anchor longest_cds -c longest_cds.conf [--overwrite] [-merge merged.fa]
                                       
    2. Using command-line arguments:
       quota_Anchor longest_cds -f sorghum.fa,maize.fa -g sorghum.gff3,maize.gff3 
                                -p sb.cds.fa,zm.cds.fa -l sorghum.cds.fa,maize.cds.fa -t 2 
                                [--overwrite] [-merge merged.fa]
                                       
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
       
       quota_Anchor longest_cds -f sorghum.fa,maize.fa -g sorghum.gff3,maize.gff3 
                                -p sb.cds.fa,zm.cds.fa -l sorghum.cds.fa,maize.cds.fa -t 2 --overwrite -c longest_cds.conf
                                [--overwrite] [-merge merged.fa]                                       
 """)
parser_sub_longest_cds.set_defaults(func=run_get_longest_cds)
parser_sub_longest_cds.add_argument('-c', '--conf', dest='conf', help="Configuration files have the lowest priority.", metavar="")
parser_sub_longest_cds.add_argument('-f', '--genome_file', dest='genome_file', help="Species genome file list(separator: ',').", metavar="")
parser_sub_longest_cds.add_argument('-g', '--gff_file', dest='gff_file', help="Species gff file list(separator: ',').", metavar="")
parser_sub_longest_cds.add_argument('-p', '--out_cds_file', dest='out_cds_file', help="Output species raw cds file list(separator: ',').", metavar="")
parser_sub_longest_cds.add_argument('-l', '--out_longest_cds_file', dest='out_longest_cds_file', help="Output species longest cds file list(separator: ',').", metavar="")
parser_sub_longest_cds.add_argument('-t', '--thread', dest='thread', help="Number of parallel processes.", type=int, metavar="")
parser_sub_longest_cds.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file.", action='store_true')
parser_sub_longest_cds.add_argument('-merge', '--merge', dest='merge', help="Optional, specify a file name, and merge all longest cds file content to this file.", metavar="")


parser_sub_get_chr_length = subparsers.add_parser('get_chr_length',
                                                  help='Generate a length file containing chromosome length and total number of genes based on the fai file and gff file.',
                                                  formatter_class=argparse.RawDescriptionHelpFormatter, description="""         
    The maize length information example file are as follows.
    chr	length	total_gene
    chr1    308452471   5892
    chr2    243675191	4751
    chr3    238017767	4103
    chr4    250330460	4093
    chr5    226353449	4485
    chr6    181357234	3412
    chr7    185808916	3070
    chr8    182411202	3536
    chr9    163004744	2988
    chr10   152435371	2705                                  
    You can execute this command in three ways: 

    1. Using configuration file:
       a)quota_Anchor get_chr_length -c [\\?, example, help] >> get_chr_length.conf 
       b)quota_Anchor get_chr_length -c get_chr_length.conf [--overwrite]    

    2. Using command-line arguments:
       quota_Anchor get_chr_length -f Oryza.fai,Sorghum.fai,Zm.fai,Setaria.fai 
            -g oryza.gff3,sorghum.gff3,maize.gff3,setaria.gff3
            -s 0-9,CHR,chr,Chr:0-9,CHR,chr,Chr:0-9,CHR,chr,Chr:0-9,CHR,chr,Chr 
            -o os_length.txt,sb_length.txt,zm_length.txt,sv_length.txt [--overwrite]    

    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.

       quota_Anchor get_chr_length -f Oryza.fai,Sorghum.fai,Zm.fai,Setaria.fai 
            -g oryza.gff3,sorghum.gff3,maize.gff3,setaria.gff3
            -s 0-9,CHR,chr,Chr:0-9,CHR,chr,Chr:0-9,CHR,chr,Chr:0-9,CHR,chr,Chr 
            -o os_length.txt,sb_length.txt,zm_length.txt,sv_length.txt 
            -c get_chr_length.conf [--overwrite]                                         
 """)
parser_sub_get_chr_length.set_defaults(func=run_lens)
parser_sub_get_chr_length.add_argument('-c', '--conf', dest='conf',
                                       help="Configuration files have the lowest priority.", metavar="")
parser_sub_get_chr_length.add_argument('-f', '--fai_file', dest='fai_file',
                                       help="Species fai files produced by gffread process or samtools(Separator: ',').",
                                       metavar="")
parser_sub_get_chr_length.add_argument('-g', '--gff_file', dest='gff_file',
                                       help="Species gff files list(Separator: ',').", metavar="")
parser_sub_get_chr_length.add_argument('-s', '--select_fai_chr_startswith', dest='select_fai_chr_startswith', help="""e.g. 0-9,chr,Chr (By default, the first column of the lines starting with numeric or chr or Chr in the fai file are extracted for plotting) (1) 0-9: software selects chromosome name start with numeric.
(2) chr: software selects chromosome name that start with 'chr'.
(3) Chr: software selects chromosome name start with 'Chr'. (Selection Separator: ',')(Species Separator: ':').""",
                                       metavar="")
parser_sub_get_chr_length.add_argument('-o', '--output_length_file', dest='output_length_file',
                                       help="Output species length file(Separator: ',').", metavar="")
parser_sub_get_chr_length.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file.",
                                       action='store_true')


# get the input file of anchorwave pro command(table file)
parser_sub_pre_col = subparsers.add_parser('pre_col', help='Generates the input file for synteny analysis (called a table file or blast file containing gene position information).', formatter_class=argparse.RawDescriptionHelpFormatter, description="""
    You can execute this command in three ways: 
                                       
    1. Using configuration file:
       a)quota_Anchor pre_col -c [\\?, example, help] >> pre_collinearity.conf 
       b)quota_Anchor pre_col -c pre_collinearity.conf [--overwrite] [--skip_blast] [-rl ref_length.txt] [-ql query_length.txt]
                                       
    2. Using command-line arguments:
       a) first method
       1)Skip blast step and use your blast file.
       quota_Anchor pre_col -b sorghum.maize.diamond -rg sorghum.gff3 -qg maize.gff3 -o sb_zm.table -bs 100 -al 0 
                                      --skip_blast [--overwrite] [-rl ref_length.txt] [-ql query_length.txt]
       2)Don't skip blast step and invoke diamond blastp.
       quota_Anchor pre_col -a diamond -rs sorghum.protein.fa -qs maize.protein.fa -db sorghum.database.diamond 
                                       -mts 20 -e 1e-10 -b sorghum.maize.diamond -rg sorghum.gff3 -qg maize.gff3
                                       -o sb_zm.table -bs 100 -al 0 
                                       [--overwrite] [-rl ref_length.txt] [-ql query_length.txt]
                                       
       b) second method
       1)Skip blast step and use your blast file.
       quota_Anchor pre_col -b sorghum.maize.blastp -rg sorghum.gff3 -qg maize.gff3 -o sb_zm.table -bs 100 -al 0 
                                      --skip_blast [--overwrite] [-rl ref_length.txt] [-ql query_length.txt]
       2)Don't skip blast step and invoke blastp.
       quota_Anchor pre_col -a blastp -rs sorghum.protein.fa -qs maize.protein.fa -db sorghum.database.blastp 
                                      -mts 20 -e 1e-10 -t 6 -ot 6 -d prot -b sorghum.maize.blastp -rg sorghum.gff3
                                      -qg maize.gff3 -o sb_zm.table -bs 100 -al 0 
                                      [--overwrite] [-rl ref_length.txt] [-ql query_length.txt]
                                      
       c) third method
       1)Skip blast step and use your blast file.
       quota_Anchor pre_col -b sorghum.maize.blastn -rg sorghum.gff3 -qg maize.gff3 -o sb_zm.table -bs 100 -al 0 
                                      --skip_blast [--overwrite] [-rl ref_length.txt] [-ql query_length.txt]
       2)Don't skip blast step and invoke blastn.
       quota_Anchor pre_col -a blastn -rs sorghum.cds.fa -qs maize.cds.fa -db sorghum.database.blastn 
                                      -mts 20 -e 1e-10 -t 6 -ot 6 -d nucl 
                                      -b sorghum.maize.blastn -rg sorghum.gff3 -qg maize.gff3 -o sb_zm.table 
                                      -bs 100 -al 0 [--overwrite] [-rl ref_length.txt] [-ql query_length.txt]

    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
       
       quota_Anchor pre_col -c pre_collinearity.conf -b sorghum.maize.diamond -rg sorghum.gff3 -qg maize.gff3 
                            -o sb_zm.table -bs 100 -al 0 --skip_blast
                            [--overwrite] [-rl ref_length.txt] [-ql query_length.txt]
       quota_Anchor pre_col -c pre_collinearity.conf -a diamond -rs sorghum.protein.fa -qs maize.protein.fa 
                            -db sorghum.database.diamond -mts 20 -e 1e-10 -b sorghum.maize.diamond 
                            -rg sorghum.gff3 -qg maize.gff3 -bs 100 -al 0 -o sb_zm.table 
                            [--overwrite] [-rl ref_length.txt] [-ql query_length.txt]                                  
 """)
parser_sub_pre_col.set_defaults(func=run_pre_col)
parser_sub_pre_col.add_argument('-c', '--conf', dest='conf', help="Configuration files have the lowest priority.", metavar="")
parser_sub_pre_col.add_argument('-skip_blast', '--skip_blast', dest='skip_blast', help="Don't run blast step and use your blast file.", action='store_true')

parser_sub_pre_col.add_argument('-a', '--align', dest='align', help="Invoking blastp/diamond/blastn alignment. If you set '-skip_blast', this parameter is useless.", choices=['diamond', 'blastp', 'blastn'], metavar="")
parser_sub_pre_col.add_argument('-rs', '--ref_seq', dest='ref_seq', help="Reference protein/cds sequence for blast. If you set '-skip_blast', this parameter is useless.", metavar="")
parser_sub_pre_col.add_argument('-qs', '--query_seq', dest='query_seq', help="Query protein/cds sequence for blast. If you set '-skip_blast', this parameter is useless.", metavar="")
parser_sub_pre_col.add_argument('-db', '--database_name', dest='database_name', help="Database name which is constructed from reference protein/cds sequence. If you set '-skip_blast', this parameter is disabled.", metavar="")
parser_sub_pre_col.add_argument('-mts', '--max_target_seqs', dest='max_target_seqs', help="Maximum number of aligned sequences to keep. If you set '-skip_blast', this parameter is disabled.", metavar="", type=int)
parser_sub_pre_col.add_argument('-e', '--evalue', dest='evalue', help="Blast e_value threshold. If you set '-skip_blast', this parameter is disabled.", metavar="")

parser_sub_pre_col.add_argument('-t', '--thread', dest='thread', help="Thread number for '-a/--align blastp/blastn'. If you set '-skip_blast', this parameter is disabled.", metavar="", type=int)
parser_sub_pre_col.add_argument('-ot', '--outfmt', dest='outfmt', help="Blast result file format(default: 6) for '-a/--align blastp/blastn'. If you set '-skip_blast', this parameter is disabled.", type=int, metavar="")
parser_sub_pre_col.add_argument('-d', '--dtype', dest='dtype', help="Sequence data type(nucl/prot) for '-a/--align blastp/blastn'. If you set '-skip_blast', this parameter is disabled.", choices=['nucl', 'prot'], metavar="")

parser_sub_pre_col.add_argument('-b', '--blast_file', dest='blast_file', help="Blast file which will be generated or had been generated.", metavar="")
parser_sub_pre_col.add_argument('-rg', '--ref_gff_file', dest='ref_gff_file', help="Reference species gff file.", metavar="")
parser_sub_pre_col.add_argument('-qg', '--query_gff_file', dest='query_gff_file', help="Query species gff file.", metavar="")
parser_sub_pre_col.add_argument('-o', '--output_file', dest='output_file', help="Output table file and this file can be used as synteny input.", metavar="")
parser_sub_pre_col.add_argument('-bs', '--bitscore', dest='bitscore', help="Filter BLAST matches to retain only those with a bitscore value greater than the specified threshold(default: 100).", metavar="", type=int)
parser_sub_pre_col.add_argument('-al', '--align_length', dest='align_length', help="Filter BLAST matches to retain only those with a alignment length value greater than the specified threshold(default: 0).", metavar="", type=int)
parser_sub_pre_col.add_argument('-rl', '--ref_length', dest='ref_length', help="Reference species length file(optional).", metavar="")
parser_sub_pre_col.add_argument('-ql', '--query_length', dest='query_length', help="Query species length file(optional).", metavar="")

parser_sub_pre_col.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file", action='store_true')


# produce collinearity file
parser_sub_col = subparsers.add_parser('col', help='Generate a collinearity file based on the table file.', formatter_class=argparse.RawDescriptionHelpFormatter, description="""
    The r_value means the maximum times of occurrences of a gene in the collinearity file, and q_value means the same as r_value.
    For maize and sorghum, maize has undergone an additional whole-genome duplication compared to sorghum.
    If sorghum is used as a reference, you can set r_value to 2 and q_value to 1.
    You can execute this command in three ways: 
                                       
    1. Using configuration file:
       a)quota_Anchor col -c [\\?, example, help] >> collinearity.conf 
       b)quota_Anchor col -c collinearity.conf [--overwrite]
                                       
    2. Using command-line arguments:
       a) get collinearity result and specify -r -q parameter
       quota_Anchor col -i sb_zm.table -o sb_zm.collinearity -r 2 -q 1 -s 0 -a 0 [--overwrite]
       
       b) get all collinearity result and remove relative inversion gene pair
       quota_Anchor col -i sb_zm.table -o sb_zm.collinearity -s 1 -a 1 [--overwrite]
       
       c) get all collinearity result and retain relative inversion gene pair
       quota_Anchor col -i sb_zm.table -o sb_zm.collinearity -s 0 -a 1 [--overwrite]    
                                         
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
                                       
       quota_Anchor col -c collinearity.conf -i sb_zm.table -o sb_zm.collinearity -r 2 -q 1 -s 1 -a 0 [--overwrite]                                      
 """)
parser_sub_col.set_defaults(func=run_col)
parser_sub_col.add_argument('-c', '--conf', dest='conf', help="Configuration files have the lowest priority.", metavar="")
parser_sub_col.add_argument('-i', '--input_file_name', dest='input_file_name', help="Input file(table file).", metavar="")
parser_sub_col.add_argument('-o', '--output_coll_name', dest='output_coll_name', help="Output collinearity file.", metavar="")
parser_sub_col.add_argument('-r', '--r_value', dest='r_value', help="Reference genome maximum alignment coverage.", metavar="", type=int)
parser_sub_col.add_argument('-q', '--q_value', dest='q_value', help="Query genome maximum alignment coverage.", metavar="", type=int)
parser_sub_col.add_argument('-s', '--strict_strand', dest='strict_strand', help="Specify whether the direction of the gene pairs within a block must be strictly the same or reverse as the block's direction(1:yes;0:no. default: 1).", metavar="", type=int)
parser_sub_col.add_argument('-a', '--get_all_collinearity', dest='get_all_collinearity', help=
                         """Enable this flag to disable r and q parameters and get all collinear result(default: 0). 
                         Options: 0: enable -r -q parameter; 1 or other integer: disable -r -q parameter and get all collinear result.""", metavar="", type=int)
parser_sub_col.add_argument('-m', '--tandem_length', dest='tandem_length',metavar="", type=int,
                            help=
                         """ This parameter is useful only for self vs self synteny alignment. Options: 0 means retain tandem gene pairs;
                             1 or any other integer means remove gene pairs with a tandem length shorter than the specified integer value(default: 0).
                             When you are doing positioning wgd events relative to species divergent events, you need set this parameter(e.g. -m 500).""")
parser_sub_col.add_argument('-I', '--minimum_chain_score', dest='minimum_chain_score', type=float, help="minimum chain score (default: 3).", metavar="")
parser_sub_col.add_argument('-W', '--overlap_window', dest='overlap_window', type=int, help="Collapse BLAST matches. Specify the maximum distance allowed, and only retain best homology pair to synteny analysis under this distance condition(default: 1).", metavar="")
parser_sub_col.add_argument('-D', '--maximum_gap_size', dest='maximum_gap_size', type=int, help="Maximum gap size for chain (default: 25).", metavar="")
parser_sub_col.add_argument('-E', '--gap_extend_penalty', dest='gap_extend_penalty', type=float, help="Chain gap extend penalty (default: -0.005).", metavar="")
parser_sub_col.add_argument('-f', '--strict_remove_overlap', dest='strict_remove_overlap', type=int, help="Specify whether to strictly remove square region gene pairs for a block to avoid overlap. (1:yes;0:no. default:  0).", metavar="")
parser_sub_col.add_argument('-t', '--count_style', dest='count_style', type=int, help="-r -q parameter's count style for a block, 0: count only the syntenic genes within a block; 1 or other integer: count all genes within a block(default: 0).", metavar="")
parser_sub_col.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file.", action='store_true')

parser_sub_ks = subparsers.add_parser('ks',
                                      help='Synonymous/non-synonymous substitution rates for syntenic gene pairs calculated in parallel.',
                                      formatter_class=argparse.RawDescriptionHelpFormatter, description="""                                           
    You can execute this command in three ways: 

    1. Using configuration file:
       a)quota_Anchor ks -c [\\?, example, help] >> ks.conf 
       b)quota_Anchor ks -c ks.conf [--overwrite] 

    2. Using command-line arguments:
       quota_Anchor ks -a mafft -i sb_zm.collinearity -p sb_zm.pep.fa
                                           -d sb_zm.cds.fa -o sb_zm.ks -t 6 [--overwrite] [--add_ks]
       quota_Anchor ks -a muscle -i zm_zm.collinearity -p zm.pep.fa
                                           -d zm.cds.fa -o zm_zm.ks -t 6[--overwrite] [--add_ks]                  
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.

       quota_Anchor ks -c ks.conf -a mafft -i sb_zm.collinearity -p sb_zm.pep.fa
                                           -d sb_zm.cds.fa -o sb_zm.ks -t 6[--overwrite] [--add_ks]
 """)
parser_sub_ks.set_defaults(func=run_ks)
parser_sub_ks.add_argument('-c', '--conf', dest='conf', help="Configuration files have the lowest priority.",
                           metavar="")
parser_sub_ks.add_argument('-i', '--collinearity', dest='collinearity', help="Collinearity file.", metavar="")
parser_sub_ks.add_argument('-a', '--align_software', dest='align_software',
                           help="Align software for every syntenic gene pair(muscle/mafft).",
                           choices=["mafft", "muscle"], metavar="")
parser_sub_ks.add_argument('-p', '--pep_file', dest='pep_file', help="Species longest protein sequence file(Separator: ',').",
                           metavar="")
parser_sub_ks.add_argument('-d', '--cds_file', dest='cds_file', help="Species longest cds sequence file(Separator: ',').", metavar="")
parser_sub_ks.add_argument('-o', '--ks_file', dest='ks_file', help="Output ks file.", metavar="")
parser_sub_ks.add_argument('-t', '--process', dest='process', help="Number of parallel processes.", metavar="",
                           type=int)
parser_sub_ks.add_argument('-add_ks', '--add_ks', dest='add_ks', help="Add extra syntenic pairs rather than overwrite it",
                           action='store_true')
parser_sub_ks.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file.",
                           action='store_true')


# collinearity dotplot or blast dotplot
parser_sub_dotplot = subparsers.add_parser('dotplot', help='Generate collinear gene pairs dotplot or homologous gene pairs dotplot.', formatter_class=argparse.RawDescriptionHelpFormatter, description="""                                           
    You can execute this command in three ways: 
                                       
    1. Using configuration file:
       a)quota_Anchor dotplot -c [\\?, example, help] >> dotplot.conf 
       b)quota_Anchor dotplot -c dotplot.conf [--overwrite] [-rm "chr"]    
                                       
    2. Using command-line arguments:
       quota_Anchor dotplot -i sb_zm.table -o sb_zm.table.png -r sb_length.txt -q zm_length.txt 
                                           -t order -r_label "Sorghum bicolor" -q_label "Zea mays"  -w 1500 -e 2000 
                                           [--overwrite] [-disable] [-use_identity] [-ks zm_sb.ks] [-a "0,3"] [-rm "chr"]   
                                       
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
      
       quota_Anchor dotplot -c dotplot.conf -i sb_zm.table -o sb_zm.table.png -r sb_length.txt -q zm_length.txt 
                                            -t order -r_label "Sorghum bicolor" -q_label "Zea mays"  -w 1500 -e 2000 
                                            [--overwrite] [-disable] [-use_identity] [-ks zm_sb.ks] [-a "0,3"] [-r "chr"]                                        
 """)
parser_sub_dotplot.set_defaults(func=run_dotplot)
parser_sub_dotplot.add_argument('-c', '--conf', dest='conf', help="Configuration files have the lowest priority.", metavar="")
parser_sub_dotplot.add_argument('-i', '--input_file', dest='input_file', help="Table file or collinearity file.", metavar="")
parser_sub_dotplot.add_argument('-o', '--output_file_name', dest='output_file_name', help="Specify a file name to save figure.", metavar="")
parser_sub_dotplot.add_argument('-r', '--ref_length', dest='ref_length', help="Reference species length file.", metavar="")
parser_sub_dotplot.add_argument('-q', '--query_length', dest='query_length', help="Query species length file.", metavar="")
parser_sub_dotplot.add_argument('-t', '--type', dest='type', help="Use gene count position within chromosome for plot(type: order) or base pair(physical distance) position within chromosome for plot(type: base)(defaults: order).", metavar="", choices=['order', 'base'])
parser_sub_dotplot.add_argument('-r_label', '--ref_name', dest='ref_name', help="Reference species coordinate axis label.", metavar="")
parser_sub_dotplot.add_argument('-q_label', '--query_name', dest='query_name', help="Query species coordinate axis label.", metavar="")
parser_sub_dotplot.add_argument('-w', '--plotnine_figure_width', dest='plotnine_figure_width', help="Plotnine module figure width (defaults: 1500)(unit: mm).", metavar="", type=int)
parser_sub_dotplot.add_argument('-e', '--plotnine_figure_height', dest='plotnine_figure_height', help="Plotnine module figure height (defaults: 1200)(unit: mm).", metavar="", type=int)
parser_sub_dotplot.add_argument('-rm', '--remove_chromosome_prefix', dest='remove_chromosome_prefix', help="Remove chromosome prefix to plot(e.g. chr,Chr,CHR)(Separator: ',').", metavar="")
parser_sub_dotplot.add_argument('-ks', '--ks', dest='ks', help="Collinearity gene pair ks file for collinearity plot(Optional).", metavar="")
parser_sub_dotplot.add_argument('-a', '--ks_area', dest='ks_area', help="ks area to plot if you specify -ks parameter for collinearity plot(Optional, default: 0,3).", metavar="")
parser_sub_dotplot.add_argument('-disable', '--disable_axis_text', dest='disable_axis_text', help="Optional, disable_axis_text for dotplot.", action='store_true')
parser_sub_dotplot.add_argument('-use_identity', '--use_identity', dest='use_identity', help="Optional, use identity as legend rather strand direction for table(blast) dotplot.", action='store_true')
parser_sub_dotplot.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file.", action='store_true')


# collinearity circle plot
parser_sub_circle = subparsers.add_parser('circle', help='Collinearity result visualization(circos).', formatter_class=argparse.RawDescriptionHelpFormatter, description="""                                           
    You can execute this command in three ways: 
                                       
    1. Using configuration file:
       a)quota_Anchor circle -c [\\?, example, help] >> circle.conf 
       b)quota_Anchor circle -c circle.conf [--overwrite] 
                                       
    2. Using command-line arguments:
       quota_Anchor circle -i sb_zm.collinearity -o sb_zm.circle.png -q zm_length.txt -r sb_length.txt 
                           -rn "Sorghum bicolor" -qn "Zea mays" -cf 7 -sf 7 -rm chr,CHR,Chr -fs 14,14 [--overwrite] 
                                       
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
     
       quota_Anchor circle -c circle.conf -i sb_zm.collinearity -o sb_zm.circle.png -q zm_length.txt -r sb_length.txt
                           -rn "Sorghum bicolor" -qn "Zea mays" -cf 7 -sf 7 -rm "chr,CHR,Chr" -fs 14,14 [--overwrite] 
 """)
parser_sub_circle.set_defaults(func=run_circle)
parser_sub_circle.add_argument('-c', '--conf', dest='conf', help="Configuration files have the lowest priority.", metavar="")
parser_sub_circle.add_argument('-i', '--input_file', dest='input_file', help="Collinearity file.", metavar="")
parser_sub_circle.add_argument('-o', '--output_file_name', dest='output_file_name', help="Specify a file name to save figure.", metavar="")
parser_sub_circle.add_argument('-r', '--ref_length', dest='ref_length', help="Reference species length file.", metavar="")
parser_sub_circle.add_argument('-q', '--query_length', dest='query_length', help="Query species length file.", metavar="")
parser_sub_circle.add_argument('-rn', '--ref_name', dest='ref_name', help="Reference species name which will as legend label.", metavar="")
parser_sub_circle.add_argument('-qn', '--query_name', dest='query_name', help="Query species name which will as legend label.", metavar="")
parser_sub_circle.add_argument('-rm', '--remove_chromosome_prefix', dest='remove_chromosome_prefix', help="Remove chromosome prefix to plot(e.g. chr,Chr,CHR)(Separator: ',').", metavar="")
parser_sub_circle.add_argument('-cf', '--chr_font_size', dest='chr_font_size', help="Chromosome name font size(defaults: 12).", metavar="", type=int)
parser_sub_circle.add_argument('-sf', '--species_name_font_size', dest='species_name_font_size', help="Species name font size(defaults: 12).", metavar="", type=int)
parser_sub_circle.add_argument('-fs', '--figsize', dest='figsize', help="Figure size(defaults: 14,14).", metavar="")
parser_sub_circle.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file.", action='store_true')


# collinearity line plot
parser_sub_line = subparsers.add_parser('line', help='Collinearity result visualization(line style).', formatter_class=argparse.RawDescriptionHelpFormatter, description="""                                           
    Figure from bottom to top (ref:oryza, query:sorghum; ref:sorghum, query:maize; ref:maize, query:setaria)
    length_file = os_length.txt, sb_length.txt, zm_length.txt, sv_length.txt
    species_name = oryza, sorghum, maize, setaria

    You can execute this command in three ways: 
                                       
    1. Using configuration file:
       a)quota_Anchor line -c [\\?, example, help] >> line.conf 
       b)quota_Anchor line -c line.conf [--overwrite] 
                                       
    2. Using command-line arguments:
       quota_Anchor line -i sb_zm.collinearity -o sb_zm.line.png -l sb_length.txt,zm_length.txt
                         -n "Sorghum bicolor,Zea mays" -rm chr,Chr,CHR -cf 7 -sf 7 -fs 14,14 [--overwrite] 
                                       
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
    
       quota_Anchor line -c line.conf -i sb_zm.collinearity -o sb_zm.line.png -l sb_length.txt,zm_length.txt
                         -n "Sorghum bicolor,Zea mays" -rm "" -cf 7 -sf 7 -fs 14,14 [--overwrite] 
 """)
parser_sub_line.set_defaults(func=run_line)
parser_sub_line.add_argument('-c', '--conf', dest='conf', help="Configuration files have the lowest priority.", metavar="")
parser_sub_line.add_argument('-i', '--input_file', dest='input_file', help="Collinearity file(e.g. file1,file2,file3)(Separator: ',').", metavar="")
parser_sub_line.add_argument('-o', '--output_file_name', dest='output_file_name', help="Specify a file name to save figure.", metavar="")
parser_sub_line.add_argument('-l', '--length_file', dest='length_file', help="Species length file list(e.g. file1,file2,file3)(Separator: ',').", metavar="")
parser_sub_line.add_argument('-n', '--species_name', dest='species_name', help="Species name list(e.g. name1,name2,name3)(Separator: ',').", metavar="")
parser_sub_line.add_argument('-rm', '--remove_chromosome_prefix', dest='remove_chromosome_prefix', help="Remove chromosome prefix to plot(e.g. chr,Chr,CHR)(Separator: ',').", metavar="")
parser_sub_line.add_argument('-cf', '--chr_font_size', dest='chr_font_size', help="Chromosome name font size(defaults: 7).", type=int, metavar="")
parser_sub_line.add_argument('-sf', '--species_name_font_size', dest='species_name_font_size', help="Species name font size(defaults: 7).", type=int, metavar="")
parser_sub_line.add_argument('-fs', '--figsize', dest='figsize', help="Figure size(defaults: 14,14).", metavar="")
parser_sub_line.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file.", action='store_true')


parser_sub_class_gene = subparsers.add_parser('class_gene', help='Genes or gene pairs are classified into whole genome duplication, tandem duplication, proximal duplication, transposed duplication, and dispersed duplication. For gene classification, there is also a single gene category (singleton) which has no homologous gene.',
                                      formatter_class=argparse.RawDescriptionHelpFormatter, description="""                                           
    You can execute this command in three ways: 

    1. Using configuration file:
       a)quota_Anchor class_gene -c [\\?, example, help] >> class_gene.conf 
       b)quota_Anchor class_gene -c class_gene.conf [--overwrite] 

    2. Using command-line arguments:
        quota_Anchor class_gene -b maize.maize.diamond -g maize.gff3 -q zm_zm.collinearity -qr bananaB_zm.collinearity
                                -o maize_classify_dir -p maize -s 1 -d 10 --overwrite -u
                                [--overwrite]
                               
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.

        quota_Anchor class_gene -b maize.maize.diamond -g maize.gff3 -q zm_zm.collinearity -qr bananaB_zm.collinearity
                                -o maize_classify_dir -p maize -s 1 -d 10 --overwrite -u -c class_gene.conf
                                [--overwrite]
 """)
parser_sub_class_gene.set_defaults(func=run_class_gene)
parser_sub_class_gene.add_argument('-c', '--conf', dest='conf', help="Configuration files have the lowest priority.", metavar="")
parser_sub_class_gene.add_argument('-b', '--query_blast', dest='query_blast', help="Focal species blast file.", metavar="")
parser_sub_class_gene.add_argument('-g', '--query_gff_file', dest='query_gff_file', help="Focal species gff file.", metavar="")
parser_sub_class_gene.add_argument('-q', '--query_query_collinearity', dest='query_query_collinearity', help="Focal species self vs self collinearity file.", metavar="")
parser_sub_class_gene.add_argument('-qr', '--query_ref_collinearity', dest='query_ref_collinearity', help="Collinearity file between focal species and outgroup species.", metavar="")
parser_sub_class_gene.add_argument('-o', '--out_directory', dest='out_directory', help="Output directory path.", metavar="")
parser_sub_class_gene.add_argument('-p', '--output_prefix', dest='output_prefix', help="Output file prefix.", metavar="")
parser_sub_class_gene.add_argument('-s', '--seg_anc', dest='seg_anc', help="Wgd/segmental genes are ancestral gene. (default: 1)[1:yes; 0: no].",choices=[0, 1], metavar="", type=int)

parser_sub_class_gene.add_argument('-u', '--unique', dest='unique', help="There is no intersection between different classifications.", action='store_true')
parser_sub_class_gene.add_argument('-d', '--proximal_max_distance', dest='proximal_max_distance', help="The maximum distance allowed for proximal genes(default: 10).", type=int, metavar="")
parser_sub_class_gene.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file.", action='store_true')


parser_sub_kde = subparsers.add_parser('kde', help='Focal species all syntenic pairs ks / block ks median histogram and gaussian kde curve.',
formatter_class = argparse.RawDescriptionHelpFormatter, description = """                                           
   You can execute this command in three ways: 

   1. Using configuration file:
      a)quota_Anchor kde -c [\\?, example, help] >> kde.conf 
      b)quota_Anchor kde -c kde.conf --overwrite

   2. Using command-line arguments:
        quota_Anchor kde -i zm_zm.collinearity -r maize.length.txt -q maize.length.txt -o zm.zm.kde.png 
                         -k zm.zm.ks [--overwrite]

   3. Using both a configuration file and command-line arguments:
      The configuration file has lower priority than other command-line parameters. 
      Parameters specified in the configuration file will be replaced by those provided via the command line.

        quota_Anchor kde -i zm_zm.collinearity -r maize.length.txt -q maize.length.txt -o zm.zm.kde.png 
                         -k zm.zm.ks -c kde.conf [--overwrite]
""")
parser_sub_kde.set_defaults(func=run_kde)
parser_sub_kde.add_argument('-c', '--conf', dest='conf', help="Configuration files have the lowest priority.", metavar="")
parser_sub_kde.add_argument('-k', '--ks_file', dest='ks_file', help="Ks file.", metavar="")
parser_sub_kde.add_argument('-i', '--collinearity_file', dest='collinearity_file', help="Collinearity file.", metavar="")
parser_sub_kde.add_argument('-r', '--ref_length', dest='ref_length', help="Reference species chromosome length file.", metavar="")
parser_sub_kde.add_argument('-q', '--query_length', dest='query_length', help="Query species chromosome length file.", metavar="")
parser_sub_kde.add_argument('-kr', '--ks_range', dest='ks_range', help="Ks range for plot(default:0,3).", metavar="")
parser_sub_kde.add_argument('-s', '--figsize', dest='figsize', help="Figure size(default:10,6.18).", metavar="")
parser_sub_kde.add_argument('-n', '--bins_number', dest='bins_number', help="histogram/curve bins number(default:100).", metavar="", type=int)
parser_sub_kde.add_argument('-o', '--output_file', dest='output_file', help="Output image name.", metavar="")
parser_sub_kde.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file.", action='store_true')

parser_sub_ks_fitting = subparsers.add_parser('kf', help='Ks fitting plot of the focal species whole genome duplication or ks fitting plot including the corrected ks peaks of species divergence events.',
formatter_class = argparse.RawDescriptionHelpFormatter, description = """                                           
   You can execute this command in three ways: 

   1. Using configuration file:
      a)quota_Anchor kf -c [\\?, example, help] >> ks_fitting.conf 
      b)quota_Anchor kf -c ks_fitting.conf [--disable_arrow] [--overwrite]

   2. Using command-line arguments:
      quota_Anchor kf -i zm_zm.collinearity -r maize.length.txt -q maize.length.txt 
                      -o zm.zm.png -k zm.zm.ks -components 2 -f maize -kr 0,2
                      [--correct_file outfile_divergent_peaks.csv] [--disable_arrow] [--overwrite]


   3. Using both a configuration file and command-line arguments:
      The configuration file has lower priority than other command-line parameters. 
      Parameters specified in the configuration file will be replaced by those provided via the command line.

      quota_Anchor kf -i zm_zm.collinearity -r maize.length.txt -q maize.length.txt
                      -o zm.zm.png -k zm.zm.ks -components 2 -f maize -kr 0,2 -c ks_fitting.conf
                      [--correct_file outfile_divergent_peaks.csv] [--disable_arrow] [--overwrite]
""")
parser_sub_ks_fitting.set_defaults(func=run_kf)
parser_sub_ks_fitting.add_argument('-c', '--conf', dest='conf', help="Configuration files have the lowest priority.", metavar="")
parser_sub_ks_fitting.add_argument('-k', '--ks_file', dest='ks_file', help="Focal species ks file(self vs self).", metavar="")
parser_sub_ks_fitting.add_argument('-i', '--collinearity_file', dest='collinearity_file', help="Focal species collinearity file(self vs self).", metavar="")
parser_sub_ks_fitting.add_argument('-r', '--ref_length', dest='ref_length', help="Reference species chromosome length file.", metavar="")
parser_sub_ks_fitting.add_argument('-q', '--query_length', dest='query_length', help="Query species chromosome length file.", metavar="")
parser_sub_ks_fitting.add_argument('-kr', '--ks_range', dest='ks_range', help="Ks range for plot(default:0,3).", metavar="")
parser_sub_ks_fitting.add_argument('-s', '--figsize', dest='figsize', help="Figure size(default:12,7).", metavar="")
parser_sub_ks_fitting.add_argument('-f', '--focal_species', dest='focal_species', help="Focal species name.", metavar="")
parser_sub_ks_fitting.add_argument('-m', '--method', dest='method', help="Typing best or mean. You can choose to use the best outgroup or the mean of all outgroups (default: best).", metavar="")
parser_sub_ks_fitting.add_argument('-components', '--components', dest='components', help="Number of whole genome duplication events and you can determine this number by dotplot.", metavar="", type=int)
parser_sub_ks_fitting.add_argument('-o', '--output_file', dest='output_file', help="Output figure file.", metavar="")
parser_sub_ks_fitting.add_argument('-correct_file', '--correct_file', dest='correct_file', help="Optional, correct file(generated by correct command).", metavar="")
parser_sub_ks_fitting.add_argument('-disable_arrow', '--disable_arrow', dest='disable_arrow', help="Disable arrow in the figure if you specify -correct_file parameter.", action='store_true')
parser_sub_ks_fitting.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file.", action='store_true')

parser_sub_trios = subparsers.add_parser('trios', help='Generate trios (consist of focal species, sister species, and outgroup species) and species pair files based on the binary tree in newick format.',
formatter_class = argparse.RawDescriptionHelpFormatter, description = """                                           
   
   You can execute this command in three ways: 

   1. Using configuration file:
      a)quota_Anchor trios -c [\\?, example, help] >> trios.conf 
      b)quota_Anchor trios -c trios.conf --overwrite

   2. Using command-line arguments:
      quota_Anchor trios -n "(((maize, sorghum), setaria), oryza);" -k maize -ot ortholog_trios_maize.csv
                         -op species_pairs.csv -t tree.txt [--overwrite]

   3. Using both a configuration file and command-line arguments:
      The configuration file has lower priority than other command-line parameters. 
      Parameters specified in the configuration file will be replaced by those provided via the command line.

      quota_Anchor trios -n "(((maize, sorghum), setaria), oryza);" -k maize -ot ortholog_trios_maize.csv
                         -op species_pairs.csv -t tree.txt -c trios.conf [--overwrite]
""")
parser_sub_trios.set_defaults(func=run_trios)
parser_sub_trios.add_argument('-c', '--conf', dest='conf', help="Configuration files have the lowest priority.", metavar="")
parser_sub_trios.add_argument('-n', '--nwk_tree', dest='nwk_tree', help="Binary tree in newick format", metavar="")
parser_sub_trios.add_argument('-k', '--focal_species', dest='focal_species', help="Focal species name within newick tree.", metavar="")
parser_sub_trios.add_argument('-ot', '--outfile_trios_path', dest='outfile_trios_path', help="Output trios file path.", metavar="")
parser_sub_trios.add_argument('-op', '--outfile_species_pair_file', dest='outfile_species_pair_file', help="Output species pair file.", metavar="")
parser_sub_trios.add_argument('-t', '--outfile_drawing_path', dest='outfile_drawing_path', help="Output newick tree file(default: tree_{focal_species}.txt).", metavar="")
parser_sub_trios.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file.", action='store_true')

parser_sub_correct = subparsers.add_parser('correct', help='The peak ks of species divergence events were fitted and corrected to the evolutionary rate level of the focal species.',
formatter_class = argparse.RawDescriptionHelpFormatter, description = """                                           
   You can execute this command in three ways: 

   1. Using configuration file:
      a)quota_Anchor correct -c [\\?, example, help] >> correct.conf 
      b)quota_Anchor correct -c correct.conf --overwrite

   2. Using command-line arguments:
      quota_Anchor correct -k "maize_sorghum.collinearity.ks,maize_setaria.collinearity.ks,
      sorghum_setaria.collinearity.ks,maize_oryza.collinearity.ks,sorghum_oryza.collinearity.ks,setaria_oryza.collinearity.ks"
      -col "maize_sorghum.collinearity,maize_setaria.collinearity,
      sorghum_setaria.collinearity,maize_oryza.collinearity,sorghum_oryza.collinearity,setaria_oryza.collinearity" 
      -s species_pairs.csv -t ortholog_trios_maize.csv -kr 0,1 -ot outfile_divergent_peaks.csv [--overwrite]

   3. Using both a configuration file and command-line arguments:
      The configuration file has lower priority than other command-line parameters. 
      Parameters specified in the configuration file will be replaced by those provided via the command line.
      
      The order of species pairs in the species pair file(specify by -s parameter/species_pair_file)
      must be consistent with the order of the ks file(specify by -k parameter/species_pair_ks_file)
      
      quota_Anchor correct -k "maize_sorghum.collinearity.ks,maize_setaria.collinearity.ks,
      sorghum_setaria.collinearity.ks,maize_oryza.collinearity.ks,sorghum_oryza.collinearity.ks,setaria_oryza.collinearity.ks"
      -col "maize_sorghum.collinearity,maize_setaria.collinearity,
      sorghum_setaria.collinearity,maize_oryza.collinearity,sorghum_oryza.collinearity,setaria_oryza.collinearity" 
      -s species_pairs.csv -t ortholog_trios_maize.csv -kr 0,1 -ot outfile_divergent_peaks.csv -c correct.conf [--overwrite]

""")

parser_sub_correct.set_defaults(func=run_correct)
parser_sub_correct.add_argument('-c', '--conf', dest='conf', help="Configuration files have the lowest priority.", metavar="")
parser_sub_correct.add_argument('-t', '--trios_file', dest='trios_file', help="Trios file based on binary tree in newick format.", metavar="")
parser_sub_correct.add_argument('-s', '--species_pair_file', dest='species_pair_file', help="Species pair file based on binary tree in newick format.", metavar="")
parser_sub_correct.add_argument('-col', '--species_pair_collinearity_file', dest='species_pair_collinearity_file', help="Collinearity file for each species pair(Separator: ','). You may need to use quotes around the option.", metavar="")
parser_sub_correct.add_argument('-k', '--species_pair_ks_file', dest='species_pair_ks_file', help="Ks file for each species pair(Separator: ','). You may need to use quotes around the option.", metavar="")
parser_sub_correct.add_argument('-kr', '--ks_range', dest='ks_range', help="Ks range for plot(default:0,3).", metavar="")
parser_sub_correct.add_argument('-ot', '--outfile_divergent_peaks', dest='outfile_divergent_peaks', help="The output file contains species ks divergent peak information.", metavar="")
parser_sub_correct.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file.", action='store_true')
#TODO: adding a space ? or parse sys.argv
args = parser.parse_args()

def copy_config_file():
    copy_dict = {
        "only_longest_pep": "longest_pep.conf",
        "pre_col": "pre_collinearity.conf",
        "longest_pep":  "longest_pep.conf",
        "longest_cds":  "longest_cds.conf",
        "col": "collinearity.conf",
        "get_chr_length": "get_chr_length.conf",
        "dotplot": "dotplot.conf",
        "circle": "circle.conf",
        "line": "line.conf",
        "ks": "ks.conf",
        "class_gene": "class_gene.conf",
        "kde": "kde.conf",
        "kf": "ks_fitting.conf",
        "trios": "trios.conf",
        "correct": "correct.conf"
    }
    return copy_dict

def main():
    # Namespace(analysis='col', conf=None, input_file=None, output_file=None, r_value=None, q_value=None, 
    # strict_strand=None, get_all_collinearity=None, count_style=None, 
    # tandem_length=None, over_lap_window=None, maximum_gap_size=None, func=<function run_coll at 0x7fe24640d620>)
    logger = logging.getLogger('main')
    logger.setLevel(level=logging.INFO)
    # logging.basicConfig(level=logging.INFO,
    #                    datefmt='%Y/%m/%d %H:%M:%S',
    #                    format='[%(asctime)s - %(name)s - %(levelname)s - %(lineno)d] - %(module)s - %(message)s')
    logging.basicConfig(level=logging.INFO,
                        # filename='output.log',
                       datefmt='%Y/%m/%d %H:%M:%S',
                       format='[%(asctime)s %(levelname)s] %(message)s')
    if hasattr(args, 'func'):
        print_help_condition = True
        for key, value in vars(args).items():
            if key != "func" and key != "analysis" and value is not None and value != False:
                print_help_condition = False
                break
        if args.conf in ["?", "help", "example"]:
            copy_dict = copy_config_file()
            config_file_path = os.path.join(base_dir, 'config_file', copy_dict[args.analysis])
            f = open(config_file_path)
            print(f.read())
        elif args.conf is None and print_help_condition:
            subparsers.choices[args.analysis].print_help()
        else:
            args.func(args)
    else:
        parser.print_help()
        sys.exit(0)
