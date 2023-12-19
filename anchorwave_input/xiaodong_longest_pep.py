import sys
import xiao_read_gff_get_longest
from Bio import SeqIO
from argparse import ArgumentParser
import re


def processing_pep(raw_pep_fil, gff_fil, longest_pep_fil):

    trans_list = xiao_read_gff_get_longest.read_gff(gff_fil)
    seqs = []
    # n = 0
    for seq_record in SeqIO.parse(raw_pep_fil, "fasta"):
        seq_record.id = seq_record.description.split(' ')[4].split(':')[1]
        seq_record.id = re.split("\\.", seq_record.id)[0]
        if seq_record.id in trans_list:
            seq_record.id = seq_record.description.split(' ')[3].split(':')[1]
            seq_record.id = re.split("\\.", seq_record.id)[0]
            # n += 1
        else:
            continue
        seqs.append(seq_record)
    SeqIO.write(seqs, longest_pep_fil, "fasta")
    # print(n)


if __name__ == '__main__':
    parser = ArgumentParser(description="This is a longest_pep from raw_pep_fil")
    parser.add_argument("-r", "--raw",
                        help="raw_pep_fil",
                        dest="raw_pep",
                        default="")
    parser.add_argument("-g", "--gff",
                        help="gff_fil",
                        dest="gff",
                        default="")
    parser.add_argument("-l", "--longest",
                        help="longest_pep_fil",
                        dest="longest",
                        default="")
    args = parser.parse_args()
    if args.raw_pep == "":
        print("Error occurs,input string is blank", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    if args.gff == "":
        print("Error occurs,please specify --ref", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    if args.longest == "":
        print("Error occurs,please specify --blast_result", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    raw_pep_file = args.raw_pep
    gff_file = args.gff
    longest_pep_file = args.longest
    processing_pep(raw_pep_file, gff_file, longest_pep_file)
