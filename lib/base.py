import pandas as pd


def read_collinearity(collinearity_file):
    coll_df = pd.read_csv(collinearity_file, header=0, comment="#", sep=",")
    # query_start = coll_df.loc[:, "queryStart"]/1000000
    # ref_start = coll_df.loc[:, "referenceStart"]/1000000
    query_id = coll_df.loc[:, "queryId"]
    ref_id = coll_df.loc[:, "refId"]
    query_chr = set()
    ref_chr = set()
    strand = coll_df.loc[:, "strand"]
    with open(collinearity_file) as f:
        _ = next(f)
        for line in f:
            if line.startswith('#'):
                query_chr.add(line.split()[4].split(sep="&")[1])
                ref_chr.add(line.split()[4].split(sep="&")[0])
    # return query_start, ref_start, query_id, ref_id, query_chr, ref_chr, strand
    return coll_df
    # collinearity_df = pd.read_table(collinearity_result, comment="#", header=0)
    #
    #
    # ref_gene_list = list(collinearity_df.loc[:, "refGene"])   # zea mays (ref) collinearity gene list
    # query_gene_list = list(collinearity_df.loc[:, "queryGene"])  # sorghum (query) collinearity gene list
    # assert (len(ref_gene_list) == len(query_gene_list))
    # for i in range(len(query_gene_list)):
    #     query_gene_list[i] = query_gene_list[i].split(sep="-")[1]
    # ref_to_query = dict(zip(ref_gene_list, query_gene_list))
    # query_to_ref = list(zip(query_gene_list, ref_gene_list))
    # return query_to_ref, block_len
    # return coll_df


def read_blast_table(combine_strand_table):
    strand_table = pd.read_csv(combine_strand_table, header=None, comment="#")
    return strand_table
