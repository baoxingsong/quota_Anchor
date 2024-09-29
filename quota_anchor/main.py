import configparser
import argparse
from pathlib import Path
import os
import sys
from .lib import base, collinearity, dotplot, circle, get_chr_length, line, line_proali_pangenome
from .lib import get_longest_pep, pre_collinearity, get_longest_cds

base_dir = Path(__file__).resolve().parent

# TODO
# software path


def run_get_longest_pep(parameter):
    global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        config_par.read(parameter.conf)
        base.file_empty(parameter.conf)
    config_soft = configparser.ConfigParser()
    config_soft.read(os.path.join(base_dir, 'config_file/software_path.ini'))
    get_longest_pep.Longest(config_par,config_soft, parameter).run_all_process()

def run_get_longest_cds(parameter):
    global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        config_par.read(parameter.conf)
        base.file_empty(parameter.conf)
    config_soft = configparser.ConfigParser()
    config_soft.read(os.path.join(base_dir, 'config_file/software_path.ini'))
    get_longest_cds.Longest(config_par,config_soft, parameter).run_all_process()

def run_pre_col(parameter):
    global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        config_par.read(parameter.conf)
        base.file_empty(parameter.conf)
    config_soft = configparser.ConfigParser()
    config_soft.read(os.path.join(base_dir, 'config_file/software_path.ini'))
    pre_collinearity.Prepare(config_par,config_soft, parameter).run_all_process()


def run_col(parameter):
    global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        config_par.read(parameter.conf)
        base.file_empty(parameter.conf)
     
    config_soft = configparser.ConfigParser()
    config_soft.read(os.path.join(base_dir, 'config_file/software_path.ini'))
    collinearity.Collinearity(config_par, config_soft, parameter).run()


def run_dotplot(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        config_par.read(parameter.conf)
        base.file_empty(parameter.conf)
    dotplot.Dotplot(config_par, parameter).run()


def run_circle(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        config_par.read(parameter.conf)
        base.file_empty(parameter.conf)
    circle.Circle(config_par, parameter).run()


def run_lens(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        config_par.read(parameter.conf)
        base.file_empty(parameter.conf)
    get_chr_length.Lens(config_par, parameter).run()


def run_line(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        config_par.read(parameter.conf)
        base.file_empty(parameter.conf)
    line.Line(config_par, parameter).run()


def run_line_3(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    base.file_empty(parameter.conf)
    line_proali_pangenome.Line(config_par).run()


#def run_class_gene(parameter):
#    # global base_dir
#    config_par = configparser.ConfigParser()
#    config_par.read(parameter.conf)
#    base.file_empty(parameter.conf)
#    if int(config_par["classification"]["type"]) == 1:
#        classification_gene.ClassGeneUnique(config_par).run()
#    else:
#        classification_gene.ClassGene(config_par).run()
#
#
#def run_clv(parameter):
#    # global base_dir
#    config_par = configparser.ConfigParser()
#    config_par.read(parameter.conf)
#    base.file_empty(parameter.conf)
#    number_gn_visualization.ClsVis(config_par).run()


parser = argparse.ArgumentParser(description='Conduct strand and WGD aware syntenic gene identification for a pair of genomes using the longest path algorithm implemented in AnchorWave.', prog="quota_Anchor")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1.1a0')


subparsers = parser.add_subparsers(title='Gene collinearity analysis', dest='analysis')

# get the longest protein sequence file
parser_sub_longest_pep = subparsers.add_parser('longest_pep', help='Get longest protein sequence file from gffread result', formatter_class=argparse.RawDescriptionHelpFormatter, description="""
    You can execute this command in three ways:
                                       
    1. Using configuration file:
       a)quota_Anchor longest_pep -c [\\?, example, help] >> longest_pep.conf 
       b)quota_Anchor longest_pep -c longest_pep.conf [--overwrite] [-merge merged.fa]
                                       
    2. Using command-line arguments:
       quota_Anchor longest_pep -f Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa,Zm-B73-REFERENCE-NAM-5.0.fa 
                                            -g Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3,Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 
                                            -p sb.p.fa,zm.p.fa -l sorghum.protein.fa,maize.protein.fa -t 2 -s [--overwrite] [-merge merged.fa]
                                       
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
       
       quota_Anchor longest_pep -f Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa,Zm-B73-REFERENCE-NAM-5.0.fa 
                                            -g Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3,Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 
                                            -p sb.p.fa,zm.p.fa -l sorghum.protein.fa,maize.protein.fa -t 2 -s  -c longest_pep.conf
                                            [--overwrite] [-merge merged.fa]                                      
 """)

parser_sub_longest_pep.set_defaults(func=run_get_longest_pep)
parser_sub_longest_pep.add_argument('-c', '--conf', dest='conf', help="Configure file which has minimum prioriy", metavar="")
parser_sub_longest_pep.add_argument('-s', '--use_s_parameter', dest='use_s_parameter', help="use '*' instead of '.' as stop codon translation in the gffread translation process. You need to set in general.", action='store_true')
parser_sub_longest_pep.add_argument('-f', '--genome_file', dest='genome_file', help="Species genome file list(separator: ',')")
parser_sub_longest_pep.add_argument('-g', '--gff_file', dest='gff_file', help="Species gff file list(separator: ',')")
parser_sub_longest_pep.add_argument('-p', '--out_pep_file', dest='out_pep_file', help="Output species raw protein file list(separator: ',')")
parser_sub_longest_pep.add_argument('-l', '--out_longest_pep_file', dest='out_longest_pep_file', help="Output species longest protein file list(separator: ',')")
parser_sub_longest_pep.add_argument('-t', '--thread', dest='thread', help="process number", type=int)
parser_sub_longest_pep.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file", action='store_true')
parser_sub_longest_pep.add_argument('-merge', '--merge', dest='merge', help="Optional, specify a file name, and merge all longest pep file content to this file")

# get the longest cds sequence file
parser_sub_longest_cds = subparsers.add_parser('longest_cds', help='Get longest cds sequence file from gffread result', formatter_class=argparse.RawDescriptionHelpFormatter, description="""
    You can execute this command in three ways: 
                                       
    1. Using configuration file:
       a)quota_Anchor longest_cds -c [\\?, example, help] >> longest_cds.conf 
       b)quota_Anchor longest_cds -c longest_cds.conf [--overwrite] [-merge merged.fa]
                                       
    2. Using command-line arguments:
       quota_Anchor longest_cds -f Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa,Zm-B73-REFERENCE-NAM-5.0.fa 
                                            -g Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3,Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 
                                            -p sb.cds.fa,zm.cds.fa -l sorghum.cds.fa,maize.cds.fa -t 6 [--overwrite] [-merge merged.fa]
                                       
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
       
       quota_Anchor longest_cds -f Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa,Zm-B73-REFERENCE-NAM-5.0.fa 
                                            -g Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3,Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 
                                            -p sb.cds.fa,zm.cds.fa -l sorghum.cds.fa,maize.cds.fa -t 6 --overwrite -c longest_cds.conf
                                            [--overwrite] [-merge merged.fa]                                       
 """)

parser_sub_longest_cds.set_defaults(func=run_get_longest_cds)
parser_sub_longest_cds.add_argument('-c', '--conf', dest='conf', help="Configure file which has minimum prioriy", metavar="")
parser_sub_longest_cds.add_argument('-f', '--genome_file', dest='genome_file', help="Species genome file list(separator: ',')")
parser_sub_longest_cds.add_argument('-g', '--gff_file', dest='gff_file', help="Species gff file list(separator: ',')")
parser_sub_longest_cds.add_argument('-p', '--out_cds_file', dest='out_cds_file', help="Output species raw cds file list(separator: ',')")
parser_sub_longest_cds.add_argument('-l', '--out_longest_cds_file', dest='out_longest_cds_file', help="Output species longest cds file list(separator: ',')")
parser_sub_longest_cds.add_argument('-t', '--thread', dest='thread', help="process number", type=int)
parser_sub_longest_cds.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file", action='store_true')
parser_sub_longest_cds.add_argument('-merge', '--merge', dest='merge', help="Optional, specify a file name, and merge all longest cds file to this file")

# get the input file of anchorwave pro command(table file)
parser_sub_pre_col = subparsers.add_parser('pre_col', help='Get input file of synteny analysis(table file)', formatter_class=argparse.RawDescriptionHelpFormatter, description="""
    You can execute this command in three ways: 
                                       
    1. Using configuration file:
       a)quota_Anchor pre_col -c [\\?, example, help] >> pre_collinearity.conf 
       b)quota_Anchor pre_col -c pre_collinearity.conf [--overwrite] [--skip_blast]
                                       
    2. Using command-line arguments:
       1)Skip blast step and use your blast file.
       quota_Anchor pre_col -b sorghum.maize.diamond -rg Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3 -qg Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -o sb_zm.table -bs 100 -al 0 
                                        --skip_blast [--overwrite]
       2)Don't skip blast step and invoke diamond blastp.
       quota_Anchor pre_col -a diamond -rs sorghum.protein.fa -qs maize.protein.fa -db sorghum.database.diamond -ob sorghum.maize.diamond -mts 20 -e 1e-10
                                       -b sorghum.maize.diamond -rg Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3 -qg Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -o sb_zm.table 
                                       -bs 100 -al 0 [--overwrite]

       1)Skip blast step and use your blast file.
       quota_Anchor pre_col -b sorghum.maize.blastp -rg Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3 -qg Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -o sb_zm.table -bs 100 -al 0 
                                        --skip_blast [--overwrite] 
       2)Don't skip blast step and invoke blastp.
       quota_Anchor pre_col -a blastp -rs sorghum.protein.fa -qs maize.protein.fa -db sorghum.database.blastp -ob sorghum.maize.blastp -mts 20 -e 1e-10 -t 6 -ot 6 -d prot 
                                      -b sorghum.maize.blastp -rg Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3 -qg Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -o sb_zm.table
                                      -bs 100 -al 0 [--overwrite]
       
       1)Skip blast step and use your blast file.
       quota_Anchor pre_col -b sorghum.maize.blastn -rg Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3 -qg Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -o sb_zm.table -bs 100 -al 0 
                                        --skip_blast [--overwrite] 
       2)Don't skip blast step and invoke blastn.
       quota_Anchor pre_col -a blastn -rs sorghum.cds.fa -qs maize.cds.fa -db sorghum.database.blastn -ob sorghum.maize.blastn -mts 20 -e 1e-10 -t 6 -ot 6 -d nucl 
                                      -b sorghum.maize.blastn -rg Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3 -qg Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -o sb_zm.table 
                                      -bs 100 -al 0 [--overwrite]

    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
       
       quota_Anchor pre_col -c pre_collinearity.conf -b sorghum.maize.diamond -rg Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3 -qg Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 
                            -o sb_zm.table -bs 100 -al 0 --skip_blast [--overwrite]
       quota_Anchor pre_col -c pre_collinearity.conf -a diamond -rs sorghum.protein.fa -qs maize.protein.fa -db sorghum.database.diamond -ob sorghum.maize.diamond -mts 20 
                            -e 1e-10 -b sorghum.maize.diamond -rg Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3 -qg Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 
                            -bs 100 -al 0 -o sb_zm.table [--overwrite]                                  
 """)

parser_sub_pre_col.set_defaults(func=run_pre_col)
parser_sub_pre_col.add_argument('-c', '--conf', dest='conf', help="Configure file which has minimum prioriy", metavar="")
parser_sub_pre_col.add_argument('-skip_blast', '--skip_blast', dest='skip_blast', help="Dont't run blast step and use your blast file", action='store_true')
parser_sub_pre_col.add_argument('-a', '--align', dest='align', help="Invoking blastp/diamond/blastn alignment. If you set '-skip_blast', this parametr is disabled.", choices=['diamond', 'blastp', 'blastn'], metavar="")
parser_sub_pre_col.add_argument('-rs', '--ref_seq', dest='ref_seq', help="Reference pretein/cds seq for blast. If you set '-skip_blast', this parametr is disabled.", metavar="")
parser_sub_pre_col.add_argument('-qs', '--query_seq', dest='query_seq', help="Query pretein/cds seq for blast. If you set '-skip_blast', this parametr is disabled.", metavar="")
parser_sub_pre_col.add_argument('-db', '--database_name', dest='database_name', help="Database name which is constructed from reference protein/cds seq. If you set '-skip_blast', this parametr is disabled.", metavar="")
parser_sub_pre_col.add_argument('-ob', '--output_blast_result', dest='output_blast_result', help="Blast result file name. If you set '-skip_blast', this parametr is disabled.", metavar="")
parser_sub_pre_col.add_argument('-mts', '--max_target_seqs', dest='max_target_seqs', help="Maximum number of aligned sequences to keep. If you set '-skip_blast', this parametr is disabled.", metavar="")
parser_sub_pre_col.add_argument('-e', '--evalue', dest='evalue', help="Blast e_value threshold. If you set '-skip_blast', this parametr is disabled.", metavar="")
parser_sub_pre_col.add_argument('-t', '--thread', dest='thread', help="Process number for '-a/--align blastp/blastn'. If you set '-skip_blast', this parametr is disabled.", metavar="", type=int)
parser_sub_pre_col.add_argument('-ot', '--outfmt', dest='outfmt', help="Blast result file format(default: 6). If you set '-skip_blast', this parametr is disabled.", type=int, metavar="")
parser_sub_pre_col.add_argument('-d', '--dtype', dest='dtype', help="Sequence data type. If you set '-skip_blast', this parametr is disabled.", choices=['nucl', 'prot'], metavar="")
parser_sub_pre_col.add_argument('-b', '--blast_file', dest='blast_file', help="Blast file which will be generated or had been generated", metavar="")
parser_sub_pre_col.add_argument('-rg', '--ref_gff_file', dest='ref_gff_file', help="Reference gff file", metavar="")
parser_sub_pre_col.add_argument('-qg', '--query_gff_file', dest='query_gff_file', help="Query gff file", metavar="")
parser_sub_pre_col.add_argument('-o', '--output_file', dest='output_file', help="Output table file and this file can be used as synteny input", metavar="")
parser_sub_pre_col.add_argument('-bs', '--bitscore', dest='bitscore', help="Filter BLAST matches to retain only those with a bitscore value greater than the specified threshold(default: 100).", metavar="", type=int)
parser_sub_pre_col.add_argument('-al', '--align_length', dest='align_length', help="Filter BLAST matches to retain only those with a alignment length value greater than the specified threshold(default: 0).", metavar="", type=int)
parser_sub_pre_col.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file", action='store_true')

# produce collinearity file
parser_sub_col = subparsers.add_parser('col', help='Get gene collinearity result file', formatter_class=argparse.RawDescriptionHelpFormatter, description="""
    You can execute this command in three ways: 
                                       
    1. Using configuration file:
       a)quota_Anchor col -c [\\?, example, help] >> collinearity.conf 
       b)quota_Anchor col -c collinearity.conf [--overwrite]
                                       
    2. Using command-line arguments:
       quota_Anchor col -i input_file_path -o output_file -r r_value -q q_value -s 1 -a 0 [--overwrite]
                                       
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
                                       
       quota_Anchor col -c collinearity.conf -i input_file_path -o output_file -r r_value -q q_value -s 1 -a 0 [--overwrite]                                      
 """)

parser_sub_col.set_defaults(func=run_col)
parser_sub_col.add_argument('-c', '--conf', dest='conf', help="Configure file which has minimum prioriy", metavar="")
parser_sub_col.add_argument('-i', '--input_file_name', dest='input_file_name', help="Input file(table file)", metavar="")
parser_sub_col.add_argument('-o', '--output_coll_name', dest='output_coll_name', help="Output collinearity file", metavar="")
parser_sub_col.add_argument('-r', '--r_value', dest='r_value', help="Reference genome maximum alignment coverage", metavar="", type=int)
parser_sub_col.add_argument('-q', '--q_value', dest='q_value', help="Query genome maximum alignment coverage", metavar="", type=int)
parser_sub_col.add_argument('-s', '--strict_strand', dest='strict_strand', help="Specify whether the direction of the gene pairs within a block must be strictly the same or reverse as the block's direction. (1:yes;0:no. default: 1)", metavar="", type=int)
parser_sub_col.add_argument('-a', '--get_all_collinearity', dest='get_all_collinearity', help=
                         """Enable this flag to get all collinear results and disable R and Q parameters (default: 0)
                            Options: 0:enable R Q parameter; 1 or other integer: get all collinear result and disable R Q parameter""", metavar="", type=int)
parser_sub_col.add_argument('-t', '--count_style', dest='count_style', type=int, help="R Q parameter's count style for a block, 0: count only the syntenic genes within a block; 1 or other integer: count all genes within a block(default: 1)", metavar="")
parser_sub_col.add_argument('-m', '--tandem_length', dest='tandem_length',metavar="", type=int,
                            help=
                         """ This parameter is useful only for self vs self synteny alignment. Options: 0 means retain tandem gene pairs;
                             1 or any other integer means remove gene pairs with a tandem length shorter than the specified integer value.(default: 0)
                             When you are doing ks peaks analysis about WGD/Divergent event, you need set this parameter""")
parser_sub_col.add_argument('-W', '--over_lap_window', dest='over_lap_window', type=int, help="Collapse BLAST matches. Specify the maximum distance allowed, and only retain best homology pair to synteny analysis under this distance condition(default: 0)", metavar="")
parser_sub_col.add_argument('-D', '--maximum_gap_size', dest='maximum_gap_size', type=int, help="Maximum gap size for chain (default: 25)", metavar="")
parser_sub_col.add_argument('-I', '--minimum_chain_score', dest='minimum_chain_score', type=int, help="minimum chain score (default: 2)", metavar="")
parser_sub_col.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file", action='store_true')


parser_sub_get_chr_length = subparsers.add_parser('get_chr_length', help='Get chromosome length and name info from fai file', formatter_class=argparse.RawDescriptionHelpFormatter, description="""                                           
    You can execute this command in three ways: 
                                       
    1. Using configuration file:
       a)quota_Anchor get_chr_length -c [\\?, example, help] >> get_chr_length.conf 
       b)quota_Anchor get_chr_length -c get_chr_length.conf [--overwrite]    
                                       
    2. Using command-line arguments:
       quota_Anchor get_chr_length -f Oryza.fai,Sorghum.fai,Zm.fai,Setari.fai 
            -g Oryza_sativa.IRGSP-1.0.59.gff3,Sorghum_bicolor.gff3,Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3,Setaria_viridis.gff3  
            -s 0-9,CHR,chr,Chr:0-9,CHR,chr,Chr:0-9,CHR,chr,Chr:0-9,CHR,chr,Chr 
            -o os_length.txt,sb_length.txt,zm_length.txt,sv_length.txt [--overwrite]    
                                       
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
                                                  
       quota_Anchor get_chr_length -f Oryza.fai,Sorghum.fai,Zm.fai,Setari.fai 
            -g Oryza_sativa.IRGSP-1.0.59.gff3,Sorghum_bicolor.gff3,Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3,Setaria_viridis.gff3  
            -s 0-9,CHR,chr,Chr:0-9,CHR,chr,Chr:0-9,CHR,chr,Chr:0-9,CHR,chr,Chr 
            -o os_length.txt,sb_length.txt,zm_length.txt,sv_length.txt 
            -c get_chr_length.conf [--overwrite]                                         
 """)
parser_sub_get_chr_length.set_defaults(func=run_lens)
parser_sub_get_chr_length.add_argument('-c', '--conf', dest='conf', help="Configure file which file has minimum prioriy", metavar="")
parser_sub_get_chr_length.add_argument('-f', '--fai_file', dest='fai_file', help="Species fai files produced by gffread(Separator: ',') ", metavar="")
parser_sub_get_chr_length.add_argument('-g', '--gff_file', dest='gff_file', help="Species gff files list(Separator: ',') ", metavar="")
parser_sub_get_chr_length.add_argument('-s', '--select_fai_chr_startswith', dest='select_fai_chr_startswith', help="""E.g. 1) 0-9: software select chromosome name start with number.
2) chr: software select chromosome name start with the string "chr".
3) Chr: software select chromosome name start with the string "Chr". (Select Separator: ',')(Species Separator: ':') """, metavar="")
parser_sub_get_chr_length.add_argument('-o', '--output_length_file', dest='output_length_file', help="Output species length file(Separator: ',')", metavar="")
parser_sub_get_chr_length.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file", action='store_true')

# collinearity dotplot or blast dotplot
parser_sub_dotplot = subparsers.add_parser('dotplot', help='Get collinearity dotplot or blast dotplot', formatter_class=argparse.RawDescriptionHelpFormatter, description="""                                           
    You can execute this command in three ways: 
                                       
    1. Using configuration file:
       a)quota_Anchor dotplot -c [\\?, example, help] >> dotplot.conf 
       b)quota_Anchor dotplot -c dotplot.conf [--overwrite]    
                                       
    2. Using command-line arguments:
       quota_Anchor dotplot -i sb_zm.table -o sb_zm.table.png -r sb_length.txt -q zm_length.txt 
                                           -t order -r_label "Sorghum bicolor" -q_label "Zea mays"  -w 1500 -e 2000 [--overwrite]    
                                       
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
      
       quota_Anchor dotplot -c dotplot.conf -i sb_zm.table -o sb_zm.table.png -r sb_length.txt -q zm_length.txt 
                                            -t order -r_label "Sorghum bicolor" -q_label "Zea mays"  -w 1500 -e 2000 [--overwrite]                                        
 """)
parser_sub_dotplot.set_defaults(func=run_dotplot)
parser_sub_dotplot.add_argument('-c', '--conf', dest='conf', help="Configure file which file has minimum prioriy", metavar="")
parser_sub_dotplot.add_argument('-i', '--input_file', dest='input_file', help="Table file or collinearity file", metavar="")
parser_sub_dotplot.add_argument('-o', '--output_file_name', dest='output_file_name', help="Specify a file name to save figure", metavar="")
parser_sub_dotplot.add_argument('-r', '--ref_length', dest='ref_length', help="Reference species length file", metavar="")
parser_sub_dotplot.add_argument('-q', '--query_length', dest='query_length', help="Query species length file", metavar="")
parser_sub_dotplot.add_argument('-t', '--type', dest='type', help="Use gene count position within chromosome for plot(type: order) or base pair(physical distance) position within chromosome for plot(type: base)(defaults: order)", metavar="", choices=['order', 'base'])   
parser_sub_dotplot.add_argument('-r_label', '--ref_name', dest='ref_name', help="Reference species coordinate axis label", metavar="")
parser_sub_dotplot.add_argument('-q_label', '--query_name', dest='query_name', help="Query species coordinate axis label", metavar="")
parser_sub_dotplot.add_argument('-w', '--plotnine_figure_width', dest='plotnine_figure_width', help="Plotnine module figure width (defaults: 1500)(unit: mm)", metavar="", type=int)
parser_sub_dotplot.add_argument('-e', '--plotnine_figure_height', dest='plotnine_figure_height', help="Plotnine module figure height (defaults: 2000)(unit: mm)", metavar="", type=int)
parser_sub_dotplot.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file", action='store_true')

# collinearity circle plot
parser_sub_circle = subparsers.add_parser('circle', help='Collinearity result visualization', formatter_class=argparse.RawDescriptionHelpFormatter, description="""                                           
    You can execute this command in three ways: 
                                       
    1. Using configuration file:
       a)quota_Anchor circle -c [\\?, example, help] >> circle.conf 
       b)quota_Anchor circle -c circle.conf [--overwrite] 
                                       
    2. Using command-line arguments:
       quota_Anchor circle -i sb_zm.collinearity -o sb_zm.circle.png -q zm_length.txt -r sb_length.txt -rn "Sorghum bicolor" -qn "Zea mays"
                                           -cf 7 -sf 7 -rm chr,CHR,Chr [--overwrite] 
                                       
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
     
       quota_Anchor circle -c circle.conf -i sb_zm.collinearity -o sb_zm.circle.png -q zm_length.txt -r sb_length.txt
                                           -rn "Sorghum bicolor" -qn "Zea mays" -cf 7 -sf 7 -rm chr,CHR,Chr [--overwrite] 
 """)
parser_sub_circle.set_defaults(func=run_circle)
parser_sub_circle.add_argument('-c', '--conf', dest='conf', help="Configure file which file has minimum prioriy", metavar="")
parser_sub_circle.add_argument('-i', '--input_file', dest='input_file', help="Collinearity file", metavar="")
parser_sub_circle.add_argument('-o', '--output_file_name', dest='output_file_name', help="Specify a file name to save figure", metavar="")
parser_sub_circle.add_argument('-r', '--ref_length', dest='ref_length', help="Reference species length file", metavar="")
parser_sub_circle.add_argument('-q', '--query_length', dest='query_length', help="Query species length file", metavar="")
parser_sub_circle.add_argument('-rn', '--ref_name', dest='ref_name', help="Reference species name")   
parser_sub_circle.add_argument('-qn', '--query_name', dest='query_name', help="Query species name")
parser_sub_circle.add_argument('-rm', '--remove_chromosome_prefix', dest='remove_chromosome_prefix', help="Remove chromosome prefix to plot(e.g. chr,Chr,CHR)(Separator: ',')")
parser_sub_circle.add_argument('-cf', '--chr_font_size', dest='chr_font_size', help="Chromosome name font size(defaults: 7)", metavar="", type=int)
parser_sub_circle.add_argument('-sf', '--species_name_font_size', dest='species_name_font_size', help="Species name font size(defaults: 7)", metavar="", type=int)
parser_sub_circle.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file", action='store_true')

# collinearity line plot
parser_sub_line_2 = subparsers.add_parser('line', help='Collinearity result visualization', formatter_class=argparse.RawDescriptionHelpFormatter, description="""                                           
    You can execute this command in three ways: 
                                       
    1. Using configuration file:
       a)quota_Anchor line -c [\\?, example, help] >> line.conf 
       b)quota_Anchor line -c line.conf [--overwrite] 
                                       
    2. Using command-line arguments:
       quota_Anchor line -i sb_zm.collinearity -o sb_zm.line.png -l sb_length.txt,zm_length.txt
                                           -n "Sorghum bicolor,Zea mays" -rm chr,Chr,CHR -cf 7 -sf 7 [--overwrite] 
                                       
    3. Using both a configuration file and command-line arguments:
       The configuration file has lower priority than other command-line parameters. 
       Parameters specified in the configuration file will be replaced by those provided via the command line.
    
       quota_Anchor line -c line.conf -i sb_zm.collinearity -o sb_zm.line.png -l sb_length.txt,zm_length.txt
                                           -n "Sorghum bicolor,Zea mays" -rm chr,Chr,CHR -cf 7 -sf 7 [--overwrite] 
 """)
parser_sub_line_2.set_defaults(func=run_line)
parser_sub_line_2.add_argument('-c', '--conf', dest='conf', help="Configure file which file has minimum prioriy", metavar="")
parser_sub_line_2.add_argument('-i', '--input_file', dest='input_file', help="Collinearity file(e.g. file1,file2,file3)(Separator: ',') ")
parser_sub_line_2.add_argument('-o', '--output_file_name', dest='output_file_name', help="Specify a file name to save figure")
parser_sub_line_2.add_argument('-l', '--length_file', dest='length_file', help="Species length file list(e.g. file1,file2,file3)(Separator: ',')")
parser_sub_line_2.add_argument('-n', '--species_name', dest='species_name', help="Species name list(e.g. name1,name2,name3)(Separator: ',')")
parser_sub_line_2.add_argument('-rm', '--remove_chromosome_prefix', dest='remove_chromosome_prefix', help="Remove chromosome prefix to plot(e.g. chr,Chr,CHR)(Separator: ',')")
parser_sub_line_2.add_argument('-cf', '--chr_font_size', dest='chr_font_size', help="Chromosome name font size(defaults: 7)", type=int)
parser_sub_line_2.add_argument('-sf', '--species_name_font_size', dest='species_name_font_size', help="Species name font size(defaults: 7)", type=int)
parser_sub_line_2.add_argument('-overwrite', '--overwrite', dest='overwrite', help="Overwrite the output file", action='store_true')

# collinearity AnchorWave proali anchors plot
parser_sub_line_proali = subparsers.add_parser('line_proali', help='Anchors file from AnchorWave proali to visualization')
parser_sub_line_proali.set_defaults(func=run_line_3)
parser_sub_line_proali.add_argument('-c', '--conf', dest='conf', help="Configure file", metavar="")


#parser_sub_clv = subparsers.add_parser('clv', help='Class gene number visualization')
#parser_sub_clv.set_defaults(func=run_clv)
#parser_sub_clv.add_argument('-c', '--conf', dest='conf', help="Command configure file", metavar="")
#parser_sub_class_gene = subparsers.add_parser('class_gene', help='Class gene as tandem, proximal, transposed, wgd/segmental, dispersed, singletons')
#parser_sub_class_gene.set_defaults(func=run_class_gene)
#parser_sub_class_gene.add_argument('-c', '--conf', dest='conf', help="Command configure file", metavar="")


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
        "line_proali": "line_proali.conf",
        "clv": "clv.conf",
        "class_gene": "class_gene.conf",
    }
    return copy_dict

def main():
    # Namespace(analysis='col', conf=None, input_file=None, output_file=None, r_value=None, q_value=None, 
    # strict_strand=None, get_all_collinearity=None, count_style=None, 
    # tandem_length=None, over_lap_window=None, maximum_gap_size=None, func=<function run_coll at 0x7fe24640d620>)
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
