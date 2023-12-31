from lib import GffFile
import sys
from argparse import ArgumentParser
# baoxing.song@pku-iaas.edu.cn


def anchorwave_quota(refGffFile, queryGffFile, blastpresult, outputFile):
    target_output = open(outputFile, 'w')
    refChromosome_gene_dict, refChromosome_gene_list, ref_GeneName_toChr_dict = GffFile.readGff(refGffFile)
    queryChromosome_gene_dict, queryChromosome_gene_list, query_GeneName_toChr_dict = GffFile.readGff(queryGffFile)

    refGeneIndex = dict()

    for refChr in refChromosome_gene_list:
        i = 0
        while i < len(refChromosome_gene_list[refChr]):
            refGeneIndex[refChromosome_gene_list[refChr][i].name] = i + 1
            i += 1

    queryGeneIndex = dict()

    for queryChr in queryChromosome_gene_list:
        i = 0
        while i < len(queryChromosome_gene_list[queryChr]):
            queryGeneIndex[queryChromosome_gene_list[queryChr][i].name] = i + 1
            i += 1

    match_pairs = set()
    with open(blastpresult) as f:
        for line in f:
            elements = line.split()
            qseqid = elements[0]
            sseqid = elements[1]
            pident = str(float(elements[2]))
            length = int(elements[3])
            mismatch = elements[4]
            gapopen = elements[5]
            qstart = elements[6]
            qend = elements[7]
            sstart = elements[8]
            send = elements[9]
            evalue = elements[10]
            bitscore = float(elements[11])
            
            match_pair = sseqid + "_" + qseqid
            if (match_pair not in match_pairs) and (bitscore > 250) and (length > 250):
                match_pairs.add(match_pair)
                target_output.write(sseqid + "\t" + ref_GeneName_toChr_dict[sseqid] + "\t" + str(refGeneIndex[sseqid]) + "\t"
                                    + str(refChromosome_gene_dict[ref_GeneName_toChr_dict[sseqid]][sseqid].start) + "\t"
                                    + str(refChromosome_gene_dict[ref_GeneName_toChr_dict[sseqid]][sseqid].end) + "\t"
                                    + refChromosome_gene_dict[ref_GeneName_toChr_dict[sseqid]][sseqid].strand + "\t"
                                    + qseqid + "\t" + query_GeneName_toChr_dict[qseqid] + "\t" + str(queryGeneIndex[qseqid]) + "\t" +
                                    str(queryChromosome_gene_dict[query_GeneName_toChr_dict[qseqid]][qseqid].start) + "\t" +
                                    str(queryChromosome_gene_dict[query_GeneName_toChr_dict[qseqid]][qseqid].end) + "\t" +
                                    queryChromosome_gene_dict[query_GeneName_toChr_dict[qseqid]][qseqid].strand + "\t" + pident + "\n")

    target_output.close()


if __name__ == '__main__':
    parser = ArgumentParser(description='Prepare file for the strand and WGD aware syntenic gene identification function implemented in AnchorWave')
    parser.add_argument("-r", "--refGFF",
                        dest="refGffFile",
                        type=str,
                        default="",
                        help="GFF file for the reference genome")
    parser.add_argument("-q", "--queryGFF",
                        dest="queryGffFile",
                        type=str,
                        default="",
                        help="GFF file for the query genome")
    parser.add_argument("-b", "--blastResult",
                        dest="blastpresult",
                        type=str,
                        default="",
                        help="BLAST result in tabular format 6")

    parser.add_argument("-o", "--output",
                        dest="outputFile",
                        type=str,
                        default="",
                        help="output file")

    args = parser.parse_args()

    if args.refGffFile == "":
        print("Error: please specify --refGFF", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.queryGffFile == "":
        print("Error: please specify --queryGFF", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.blastpresult == "":
        print("Error: please specify --blastResult", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.outputFile == "":
        print("Error: please specify --output", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    refGffFile = args.refGffFile
    queryGffFile = args.queryGffFile
    blastpresult = args.blastpresult
    outputFile = args.outputFile
    anchorwave_quota(refGffFile, queryGffFile, blastpresult, outputFile)
