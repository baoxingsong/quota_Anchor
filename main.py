import configparser
import argparse
from pathlib import Path
import os
from .lib import pre_collinearity, collinearity, dotplot, prepare_ks, ks, blockinfo, ks_peaks, peaksfit, ksfigure, \
     classification_gene, orthogroup3, number_gn_visualization, orthogroup4, duplicate_pair_ks, circle, line, get_chr_length


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


def run_prepare_ks(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    config_soft = configparser.ConfigParser()
    config_soft.read(parameter.conf)
    prepare_ks.Prepare(config_par, config_soft).run()


def run_ks(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    config_soft = configparser.ConfigParser()
    config_soft.read(os.path.join(base_dir, 'config_file/software_path.ini'))
    ks.Ks(config_par, config_soft).first_layer_run()


def run_block_info(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    blockinfo.BlockInfo(config_par).run()


def run_ks_pdf(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    ks_peaks.KsPeaks(config_par).run()


def run_pdf_fit(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    peaksfit.PeaksFit(config_par).run()


def run_ks_figure(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    ksfigure.KsFigure(config_par).run()


# def run_pre_classification():
#     config_par = configparser.ConfigParser()
#     config_par.read('./config_file/classification.conf')
#     pre_classification.PreClassification(config_par).run()


# def run_classification():
#     config_par = configparser.ConfigParser()
#     config_par.read('./config_file/classification.conf')
#     classification.Classification(config_par).run()


def run_class_gene(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    if int(config_par["classification"]["type"]) == 1:
        classification_gene.ClassGeneUnique(config_par).run()
    else:
        classification_gene.ClassGene(config_par).run()


# def run_group():
#     global base_dir
#     config_par = configparser.ConfigParser()
#     config_par.read(os.path.join(base_dir, 'config_file/orthogroup.conf'))
#     orthogroup.Group(config_par).run()


# def run_group2():
#     global base_dir
#     config_par = configparser.ConfigParser()
#     config_par.read(os.path.join(base_dir, 'config_file/orthogroup2.conf'))
#     orthogroup2.Group(config_par).run()


def run_group3(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    orthogroup3.Group(config_par).run()


def run_group4(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    orthogroup4.Group(config_par).run()


def run_clv(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    number_gn_visualization.ClsVis(config_par).run()


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


def run_line(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    line.Line(config_par).run()


def run_dup(parameter):
    # global base_dir
    config_par = configparser.ConfigParser()
    config_par.read(parameter.conf)
    config_soft = configparser.ConfigParser()
    config_soft.read(os.path.join(base_dir, 'config_file/software_path.ini'))
    duplicate_pair_ks.DuplicateKs(config_par, config_soft).run()
# def run_class_gene_2():
#     config_par = configparser.ConfigParser()
#     config_par.read('./config_file/classification_gene_2.conf')
#     if int(config_par["classification"]["type"]) == 1:
#         classification_gene_2.ClassGeneUnique(config_par).run()
#     else:
#         classification_gene_2.ClassGene(config_par).run()


parser = argparse.ArgumentParser(description='collinearity gene analysis', prog="quota_Anchor")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1')

subparsers = parser.add_subparsers(title='gene collinearity analysis', dest='analysis')
# get the longest protein and AnchorWave pro input file
parser_sub_a = subparsers.add_parser('pre_col', help='get longest protein and AnchorWave pro(collinearity) input file')
parser_sub_a.set_defaults(func=run_pre_coll)
parser_sub_a.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
# produce collinearity file
parser_sub_b = subparsers.add_parser('col', help='get gene collinearity file by AnchorWave pro command')
parser_sub_b.set_defaults(func=run_coll)
parser_sub_b.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
# collinearity dotplot or blast dotplot
parser_sub_c = subparsers.add_parser('dotplot', help='collinearity dotplot or blast dotplot')
parser_sub_c.set_defaults(func=run_dotplot)
parser_sub_c.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
# prepare data for synonymous mutation and non-synonymous mutation(longest cds)
parser_sub_d = subparsers.add_parser('pre_ks', help='get longest cds')
parser_sub_d.set_defaults(func=run_prepare_ks)
parser_sub_d.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
# synonymous mutation and non-synonymous mutation
parser_sub_e = subparsers.add_parser('ks', help='get ks and ka information')
parser_sub_e.set_defaults(func=run_ks)
parser_sub_e.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
# summary block information
parser_sub_f = subparsers.add_parser('block_info', help='summary block information')
parser_sub_f.set_defaults(func=run_block_info)
parser_sub_f.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
# ks probability density function curve
parser_sub_g = subparsers.add_parser('kp', help='ks(total, average, median) probability density function curve')
parser_sub_g.set_defaults(func=run_ks_pdf)
parser_sub_g.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
# ks probability density function curve fitting
parser_sub_h = subparsers.add_parser('pf', help='ks probability density function curve fitting')
parser_sub_h.set_defaults(func=run_pdf_fit)
parser_sub_h.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
# ks probability density function curve fitting
parser_sub_i = subparsers.add_parser('kf', help='ks distribution figure')
parser_sub_i.set_defaults(func=run_ks_figure)
parser_sub_i.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
# #
# parser_sub1 = subparsers1.add_parser('get_chr_info', help='assumed ancestor file')
# parser_sub1.set_defaults(func=run_pre_classification)
# # get ancestor chr 1 chr_total_gene_number class
# parser_sub1 = subparsers1.add_parser('class', help='relative classification')
# parser_sub1.set_defaults(func=run_classification)
# class gene
parser_sub_j = subparsers.add_parser('class_gene', help='class gene as tandem, proximal, transposed, wgd/segmental, dispersed, singletons')
parser_sub_j.set_defaults(func=run_class_gene)
parser_sub_j.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
# parser_sub1 = subparsers1.add_parser('group', help='orthogroup based collinearity')
# parser_sub1.set_defaults(func=run_group)
# parser_sub1 = subparsers1.add_parser('group2', help='orthogroup based collinearity')
# parser_sub1.set_defaults(func=run_group2)
parser_sub_k = subparsers.add_parser('group3', help='orthogroup based collinearity')
parser_sub_k.set_defaults(func=run_group3)
parser_sub_k.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
parser_sub_l = subparsers.add_parser('group4', help='orthogroup based collinearity')
parser_sub_l.set_defaults(func=run_group4)
parser_sub_l.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
parser_sub_m = subparsers.add_parser('clv', help='class gene number visualization')
parser_sub_m.set_defaults(func=run_clv)
parser_sub_m.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
parser_sub_n = subparsers.add_parser('dupks', help='class gene number visualization')
parser_sub_n.set_defaults(func=run_dup)
parser_sub_n.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
parser_sub_o = subparsers.add_parser('circle', help='AnchorWave pro command collinearity visualization')
parser_sub_o.set_defaults(func=run_circle)
parser_sub_o.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
parser_sub_p_pre = subparsers.add_parser('get_chr_length', help='chromosome length and name info from fai file')
parser_sub_p_pre.set_defaults(func=run_lens)
parser_sub_p_pre.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
parser_sub_p = subparsers.add_parser('line', help='AnchorWave pro command collinearity visualization')
parser_sub_p.set_defaults(func=run_line)
parser_sub_p.add_argument('-c', '--conf', dest='conf', help="command configure file", metavar="")
#    parser_sub1 = subparsers1.add_parser('class_gene_2', help='class gene as tandem, proximal, transposed, wgd/segmental, dispersed, singletons')
#    parser_sub1.set_defaults(func=run_class_gene_2)

args = parser.parse_args()


def main():
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
