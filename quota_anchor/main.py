import configparser
import argparse
from pathlib import Path
import os
from .lib import pre_collinearity, collinearity, dotplot, circle, get_chr_length, line_2, line_proali_pangenome


base_dir = Path(__file__).resolve().parent

# TODO
# software path


def run_pre_coll(parameter):
    global base_dir
    config_par = configparser.ConfigParser()
    # config_par.read(os.path.join(base_dir, "config_file/pre_collinearity.conf"))
    config_par.read(parameter.conf)
    config_soft = configparser.ConfigParser()
    config_soft.read(os.path.join(base_dir, 'config_file/software_path.ini'))
    pre_collinearity.Prepare(config_par, config_soft).run_all_process()


def run_coll(parameter):
    global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    config_soft = configparser.ConfigParser()
    config_soft.read(os.path.join(base_dir, 'config_file/software_path.ini'))
    collinearity.Collinearity(config_par, config_soft).run()


def run_dotplot(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    dotplot.Dotplot(config_par).run()


def run_circle(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    circle.Circle(config_par).run()


def run_lens(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    get_chr_length.Lens(config_par).run()


def run_line_2(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    line_2.Line(config_par).run()


def run_line_3(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    line_proali_pangenome.Line(config_par).run()


parser = argparse.ArgumentParser(description='Conduct strand and WGD aware syntenic gene identification for a pair of genomes using the longest path algorithm implemented in AnchorWave', prog="quota_Anchor")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1a0')


subparsers = parser.add_subparsers(title='Gene collinearity analysis', dest='analysis')
# get the longest protein and AnchorWave pro input file
parser_sub_a = subparsers.add_parser('pre_col', help='Get longest protein file from gffread result and input file for gene collinearity analysis')
parser_sub_a.set_defaults(func=run_pre_coll)
parser_sub_a.add_argument('-c', '--conf', dest='conf', help="configure file", metavar="")
# produce collinearity file
parser_sub_b = subparsers.add_parser('col', help='Get gene collinearity result file')
parser_sub_b.set_defaults(func=run_coll)
parser_sub_b.add_argument('-c', '--conf', dest='conf', help="configure file", metavar="")


parser_sub_p_pre = subparsers.add_parser('get_chr_length', help='Get chromosome length and name info from fai file')
parser_sub_p_pre.set_defaults(func=run_lens)
parser_sub_p_pre.add_argument('-c', '--conf', dest='conf', help="configure file", metavar="")
# collinearity dotplot or blast dotplot
parser_sub_c = subparsers.add_parser('dotplot', help='Get collinearity dotplot or blast dotplot')
parser_sub_c.set_defaults(func=run_dotplot)
parser_sub_c.add_argument('-c', '--conf', dest='conf', help="configure file", metavar="")
# collinearity circle plot
parser_sub_o = subparsers.add_parser('circle', help='Collinearity result visualization')
parser_sub_o.set_defaults(func=run_circle)
parser_sub_o.add_argument('-c', '--conf', dest='conf', help="configure file", metavar="")
# collinearity line plot
parser_sub_p_2 = subparsers.add_parser('line_2', help='Collinearity result visualization')
parser_sub_p_2.set_defaults(func=run_line_2)
parser_sub_p_2.add_argument('-c', '--conf', dest='conf', help="configure file", metavar="")
# collinearity AnchorWave proali anchors plot
parser_sub_p_3 = subparsers.add_parser('line_proali', help='Anchors file from AnchorWave proali to visualization')
parser_sub_p_3.set_defaults(func=run_line_3)
parser_sub_p_3.add_argument('-c', '--conf', dest='conf', help="configure file", metavar="")


args = parser.parse_args()


def main():
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
