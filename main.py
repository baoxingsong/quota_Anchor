# import sys
import configparser
import argparse
from lib import collinearity, dotplot, ks, prepare_ks


def run_coll():
    config_par = configparser.ConfigParser()
    config_par.read('./config_file/collinearity.conf')
    config_soft = configparser.ConfigParser()
    config_soft.read('./config_file/software_path.ini')
    collinearity.Collinearity(config_par, config_soft).run_all_processes()


def run_blast_coll():
    config_par = configparser.ConfigParser()
    config_par.read('./config_file/collinearity.conf')
    config_soft = configparser.ConfigParser()
    config_soft.read('./config_file/software_path.ini')
    collinearity.Collinearity(config_par, config_soft).run_all_processes_blastp()


def run_dotplot():
    config_par = configparser.ConfigParser()
    config_par.read('./config_file/dotplot.conf')
    dotplot.Dotplot(config_par).run_coll_dotplot()


def run_blast_dotplot():
    config_par = configparser.ConfigParser()
    config_par.read('./config_file/dotplot.conf')
    dotplot.BlastDotplot(config_par).run_blast_dotplot()


def run_prepare_ks():
    config_par = configparser.ConfigParser()
    config_par.read('./config_file/ks.conf')
    config_soft = configparser.ConfigParser()
    config_soft.read('./config_file/software_path.ini')
    prepare_ks.Prepare(config_par, config_soft).run()


def run_ks():
    config_par = configparser.ConfigParser()
    config_par.read('./config_file/ks.conf')
    config_soft = configparser.ConfigParser()
    config_soft.read('./config_file/software_path.ini')
    ks.Ks(config_par, config_soft).sub_run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='collinearity gene analysis')

    subparsers1 = parser.add_subparsers(title='gene collinearity analysis', dest='collinearity_analysis')

    # produce collinearity file
    parser_sub1 = subparsers1.add_parser('col', help='get gene collinearity file by AnchorWave pro command')
    parser_sub1.set_defaults(func=run_coll)
    # produce collinearity file(blast makedb)
    parser_sub1 = subparsers1.add_parser('blast_col', help='blastp make database and alignment, get gene collinearity file by AnchorWave pro command')
    parser_sub1.set_defaults(func=run_blast_coll)
    # collinearity dotplot
    parser_sub1 = subparsers1.add_parser('dotplot', help='dot plot by AnchorWave pro command')
    parser_sub1.set_defaults(func=run_dotplot)
    # blast dotplot
    parser_sub1 = subparsers1.add_parser('blast_dotplot', help='blast dot plot by combineBlastAndStrandInformation script')
    parser_sub1.set_defaults(func=run_blast_dotplot)
    # prepare data for synonymous mutation and non-synonymous mutation(longest cds)
    parser_sub1 = subparsers1.add_parser('pre_ks', help='get longest cds')
    parser_sub1.set_defaults(func=run_prepare_ks)
    # synonymous mutation and non-synonymous mutation
    parser_sub1 = subparsers1.add_parser('ks', help='get ks and ka information')
    parser_sub1.set_defaults(func=run_ks)

    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func()
    else:
        parser.print_help()
