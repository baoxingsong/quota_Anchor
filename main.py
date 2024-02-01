import configparser
import argparse
from lib import pre_collinearity, collinearity, dotplot, prepare_ks, ks


def run_pre_coll():
    config_par = configparser.ConfigParser()
    config_par.read('./config_file/collinearity.conf')
    config_soft = configparser.ConfigParser()
    config_soft.read('./config_file/software_path.ini')
    pre_collinearity.Prepare(config_par, config_soft).run_all_process()


def run_coll():
    config_par = configparser.ConfigParser()
    config_par.read('./config_file/collinearity.conf')
    config_soft = configparser.ConfigParser()
    config_soft.read('./config_file/software_path.ini')
    collinearity.Collinearity(config_par, config_soft).run()


def run_dotplot():
    config_par = configparser.ConfigParser()
    config_par.read('./config_file/dotplot.conf')
    dotplot.Dotplot(config_par).run()


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

    # get longest protein and AnchorWave pro input file
    parser_sub1 = subparsers1.add_parser('pre_col', help='get longest protein and AnchorWave pro(collinearity) input file')
    parser_sub1.set_defaults(func=run_pre_coll)
    # produce collinearity file
    parser_sub1 = subparsers1.add_parser('col', help='get gene collinearity file by AnchorWave pro command')
    parser_sub1.set_defaults(func=run_coll)
    # collinearity dotplot or blast dotplot
    parser_sub1 = subparsers1.add_parser('dotplot', help='collinearity dotplot or blast dotplot')
    parser_sub1.set_defaults(func=run_dotplot)
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
