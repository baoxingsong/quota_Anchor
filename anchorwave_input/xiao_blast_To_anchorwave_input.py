import xiao_read_gff_get_longest
from argparse import ArgumentParser
import sys


def generate_anchor_wave_input(query_gff_fil, ref_gff_fil, blast_file, anchor_wave_input):
    query_dict_chr_gene_gene_class, query_gene_to_chr = xiao_read_gff_get_longest.read_gff_get_gene_(query_gff_fil)
    ref_dict_chr_gene_gene_class, ref_gene_to_chr = xiao_read_gff_get_longest.read_gff_get_gene_(ref_gff_fil)
    query_index_dict = dict()
    ref_index_dict = dict()
    for q in query_dict_chr_gene_gene_class:
        i = 0
        while i < len(query_dict_chr_gene_gene_class[q]):
            query_index_dict[query_dict_chr_gene_gene_class[q][i][0].name] = i+1
            i += 1
    for s in ref_dict_chr_gene_gene_class:
        i = 0
        while i < len(ref_dict_chr_gene_gene_class[s]):
            ref_index_dict[ref_dict_chr_gene_gene_class[s][i][0].name] = i+1
            i += 1
    out_ = open(anchor_wave_input, "w")
    with open(blast_file) as f:
        for line in f:
            elements = line.split()
            query_name = elements[0]
            ref_name = elements[1]
            identity = elements[2]
            alignment_length = int(elements[3])
            bit_score = float(elements[11])
            if (alignment_length > 250) and (bit_score > 250):
                out_.write(
                    query_name + "\t" + str(query_gene_to_chr[query_name]) + "\t" +
                    str(query_index_dict[query_name]) + "\t" +
                    str(query_dict_chr_gene_gene_class[query_gene_to_chr[query_name]][query_index_dict[query_name]-1][0].start) + "\t" +
                    str(query_dict_chr_gene_gene_class[query_gene_to_chr[query_name]][query_index_dict[query_name]-1][0].end) + "\t" +
                    str(query_dict_chr_gene_gene_class[query_gene_to_chr[query_name]][query_index_dict[query_name]-1][0].strand) + "\t" +
                    ref_name + "\t" + str(ref_gene_to_chr[ref_name]) + "\t" +
                    str(ref_index_dict[ref_name]) + "\t" +
                    str(ref_dict_chr_gene_gene_class[ref_gene_to_chr[ref_name]][ref_index_dict[ref_name]-1][0].start) + "\t" +
                    str(ref_dict_chr_gene_gene_class[ref_gene_to_chr[ref_name]][ref_index_dict[ref_name]-1][0].end) + "\t" +
                    str(ref_dict_chr_gene_gene_class[ref_gene_to_chr[ref_name]][ref_index_dict[ref_name]-1][0].strand) + "\t" +
                    str(identity) + "\n")
    out_.close()


if __name__ == '__main__':
    parser = ArgumentParser(description="This is a format change about collinearity")
    parser.add_argument("-q", "--query",
                        help="query gff_file",
                        dest="query_gff",
                        default="")
    parser.add_argument("-r", "--ref",
                        help="ref_gff_file",
                        dest="ref_gff",
                        default="")
    parser.add_argument("-b", "--blast_result",
                        help="blast_result gff_file",
                        dest="blast_result",
                        default="")
    parser.add_argument("-o", "--output",
                        help="outputfile",
                        dest="output",
                        default="")
    args = parser.parse_args()
    if args.query_gff == "":
        print("Error occurs,input string is blank", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    if args.ref_gff == "":
        print("Error occurs,please specify --ref", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    if args.blast_result == "":
        print("Error occurs,please specify --blast_result", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    if args.output == "":
        print("Error occurs,input string is blank", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    query_gff_file = args.query_gff
    ref_gff_file = args.ref_gff
    blast_results = args.blast_result
    output = args.output
    generate_anchor_wave_input(query_gff_file, ref_gff_file, blast_results, output)
