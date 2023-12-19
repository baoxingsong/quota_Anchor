import re
import numpy as np


class Transcript:
    pep_name = ""
    trans_name = ""
    start = 0
    end = 0
    distance = 0

    def __init__(self, star, en, trans_nam, distance):
        self.distance = distance
        self.start = star
        self.end = en
        self.trans_name = trans_nam

    def __lt__(self, other):
        if self.distance < other.distance:
            return True
        elif self.distance == other.distance and self.start < other.start:
            return True
        else:
            return False

    def __gt__(self, other):
        if self.distance > other.distance:
            return True
        elif self.distance == other.distance and self.start > other.start:
            return True
        else:
            return False

    def __eq__(self, other):
        if self.distance == other.distance and self.start == other.start:
            return True
        else:
            return False


class Gene:
    gene_name = ""
    start = 0
    end = 0
    ch = ""
    strand = ""

    def __init__(self, ch, start, end, strand, name):
        self.chr = ch
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name

    def __lt__(self, other):
        if self.start < other.start:
            return True
        elif self.start == other.start and self.end < other.end:
            return True
        else:
            return False

    def __gt__(self, other):
        if self.start > other.start:
            return True
        elif self.start == other.start and self.end > other.end:
            return True
        else:
            return False

    def __eq__(self, other):
        if self.start == other.start and self.end == other.end:
            return True
        else:
            return False


def read_gff(gff_file):
    dict_gene_transcript_pep = dict()
    # chr_trans_pep_map = dict()
    trans_list = []
    # notice here , gff file is input ///////////////////////////////////////////////////////////////////////////
    # notice here , gff file is input ///////////////////////////////////////////////////////////////////////////
    # notice here , gff file is input ///////////////////////////////////////////////////////////////////////////
    with open(gff_file) as f:
        for line in f:
            # dict_gene_transcript_pep is chr : gene_name : Transcript(class), Transcript has attribute start, end, distance, trans_name.
            m = re.search('^(\\S+)\t(\\S+)\tmRNA\t(\\S+)\t(\\S+)\t(\\S+)\t(\\S+)\t(\\S)\tID=transcript:(.+?);Parent=gene:(.+?);' or
                          "^(\\S+)\t(\\S+)\tV_gene_segment\t(\\S+)\t(\\S+)\t(\\S+)\t(\\S+)\t(\\S)\tID=transcript:(.+?);Parent=gene:(.+?);", line)
            if m is not None:
                ch = str(m.group(1))
                start = int(m.group(3))
                end = int(m.group(4))
                trans_name = str(m.group(8))
                gene_name = str(m.group(9))
                if ch not in dict_gene_transcript_pep:
                    dict_gene_transcript_pep[ch] = dict()
                else:
                    if gene_name not in dict_gene_transcript_pep[ch]:
                        dict_gene_transcript_pep[ch][gene_name] = np.empty([0, 1], dtype=Transcript)
                    transcript = Transcript(start, end, trans_name, end-start)
                    dict_gene_transcript_pep[ch][gene_name] = np.append(dict_gene_transcript_pep[ch][gene_name], np.array([[transcript]]), axis=0)
        for ch in dict_gene_transcript_pep:
            # print("aaaa")
            for ge in dict_gene_transcript_pep[ch]:
                # print(type(dict_gene_transcript_pep[ch][ge][0]))
                # print("bbbb")
                # i = 0
                # print(len(dict_gene_transcript_pep[ch][ge]))
                # print(len(dict_gene_transcript_pep[ch][ge]))
                dict_gene_transcript_pep[ch][ge] = np.sort(dict_gene_transcript_pep[ch][ge], axis=0)
                trans_list.append(dict_gene_transcript_pep[ch][ge][-1][0].trans_name)
                # print(trans_list[i])
                # print(len(trans_list))
                # i +=1
    return trans_list
    # return None
# else:
#     # ch_trans_pep_map is ch : trans_name : pep_name.
#     m = re.search("^(\\S+)\t(\\S+)\tCDS\t(\\S+)\t(\\S+)\t(\\S+)\t(\\S+)\t(\\S)\tID=CDS:(.+?);Parent=transcript:(.+?);", line)
#     if m is not None:
#         ch = str(m.group(1))
#         pep_name = str(m.group(8))
#         trans_name = str(m.group(9))
#         if ch not in ch_trans_pep_map:
#             ch_trans_pep_map[ch] = dict()
#         if trans_name not in ch_trans_pep_map[ch]:
#             ch_trans_pep_map[ch][trans_name] = pep_name


# the function is to get gene_Gene which includes gene's start, end and strand information.
def read_gff_get_gene_(gff_file):
    dict_chr_gene_gene_class = dict()
    gene_to_chr = dict()
    with open(gff_file) as f:
        for line in f:
            # dict_gene_transcript_pep is chr : gene_name : Transcript(class), Transcript has attribute start, end, distance, trans_name.
            m = re.search('^(\\S+)\t(\\S+)\tgene\t(\\S+)\t(\\S+)\t(\\S+)\t(\\S+)\t(\\S)\tID=gene:(.+?);', line)
            if m is not None:
                ch = str(m.group(1))
                start = int(m.group(3))
                end = int(m.group(4))
                strand = str(m.group(6))
                gene_name = str(m.group(8))
                if gene_name not in gene_to_chr:
                    gene_to_chr[gene_name] = ch
                if ch not in dict_chr_gene_gene_class:
                    dict_chr_gene_gene_class[ch] = np.empty([0, 1], dtype=Gene)
                gene = Gene(ch, start, end, strand, gene_name)
                dict_chr_gene_gene_class[ch] = np.append(dict_chr_gene_gene_class[ch], np.array([[gene]]), axis=0)
    for ch in dict_chr_gene_gene_class:
        dict_chr_gene_gene_class[ch] = np.sort(dict_chr_gene_gene_class[ch], axis=0)
    return dict_chr_gene_gene_class, gene_to_chr


# if __name__ == '__main__':
