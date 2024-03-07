import pandas as pd
import numpy as np


def read_ks(file, col):
    ks = pd.read_csv(file, sep='\t', low_memory=False)
    ks.drop_duplicates(subset=['id1', 'id2'], keep='first', inplace=True)
    ks[col] = ks[col].astype(float)
    ks = ks[ks[col] >= 0]
    ks.index = ks['id1'] + ',' + ks['id2']
    return ks[col]


def get_median(data):
    data = [i for i in data if i >= 0]
    if len(data) == 0:
        return 0
    data.sort()
    half = len(data) // 2
    return (data[half] + data[~half]) / 2


def get_average(data):
    data = [i for i in data if i >= 0]
    if len(data) == 0:
        return 0
    average = sum(data) / len(data)
    return average


def record_block_region(start1, start2, end1, end2, block_order):
    start1.append(block_order[0][0])
    start2.append(block_order[0][1])
    end1.append(block_order[-1][0])
    end2.append(block_order[-1][1])
    return start1, start2, end1, end2, block_order


def use_underscore_concatenate(block1, block2, ks, block_ks, block_order):
    query = [str(item[0]) for item in block_order]
    block1.append("_".join(query))
    ref = [str(item[1]) for item in block_order]
    block2.append("_".join(ref))
    block_ks = [str(i) for i in block_ks]
    ks.append("_".join(block_ks))
    return block1, block2, ks, block_ks, block_order


def read_table_file(table_file):
    table_df = pd.read_csv(table_file, header=None, sep="\t", low_memory=False)
    table_df.rename(columns={2: 'chr2_order', 8: 'chr1_order'}, inplace=True)
    table_df[0] = table_df[0].astype(str)
    table_df[6] = table_df[6].astype(str)
    map2 = table_df["chr2_order"].copy()
    map2.index = table_df[0]
    duplicate_index = map2.index.duplicated()
    map2 = map2[~duplicate_index]
    return table_df, map2


class BlockInfo:
    def __init__(self, config_par):
        self.score = int(config_par["block_info"]["score"])
        self.evalue = float(config_par["block_info"]["evalue"])

        self.blast_file = config_par["block_info"]["blast"]
        self.table_file = config_par["block_info"]["table_file"]
        self.repeat_number = int(config_par["block_info"]["repeat_number"])
        self.ks_col = config_par["block_info"]["ks_col"]
        self.ks = config_par["block_info"]["ks"]
        self.collinearity = config_par["block_info"]["collinearity"]
        self.savefile = config_par["block_info"]["savefile"]

    def tandem_ratio(self, table_df, map2, block):
        block['order2'] = block['id2'].map(map2)
        block_blast = table_df[(table_df[6].isin(block['id1'].values.tolist())) & (
            table_df[0].isin(block['id2'].values.tolist()))].copy()
        block_blast = pd.merge(
            block_blast, block, left_on=6, right_on='id1', how='left')
        block_blast['difference'] = (
            block_blast['chr2_order'] - block_blast['order2']).abs()
        block_blast = block_blast[(block_blast['difference'] <= self.repeat_number) & (
            block_blast['difference'] > 0)]
        return len(block_blast[6].unique())/len(block)*len(block_blast)/(len(block)+len(block_blast))

    def blast_homo(self, blast_file, repeat_number):
        blast = pd.read_csv(blast_file, sep="\t", header=None, low_memory=False)
        blast = blast[(blast[11] >= self.score) & (
                blast[10] < self.evalue) & (blast[1] != blast[0])]
        blast.drop_duplicates(subset=[0, 1], keep='first', inplace=True)
        blast[0] = blast[0].astype(str)
        blast[1] = blast[1].astype(str)
        index = [group.sort_values(by=11, ascending=False)[:repeat_number].index.tolist() for name, group in blast.groupby([0])]
        blast = blast.loc[np.concatenate(
            np.array([k[:repeat_number] for k in index], dtype=object)), [0, 1]]
        blast = blast.assign(homo1=np.nan, homo2=np.nan, homo3=np.nan, homo4=np.nan, homo5=np.nan)
        for i in range(1, 6):
            blue_num = i + 5
            red_index = np.concatenate(
                np.array([k[:i] for k in index], dtype=object))
            blue_index = np.concatenate(
                np.array([k[i:blue_num] for k in index], dtype=object))
            gray_index = np.concatenate(
                np.array([k[blue_num:repeat_number] for k in index], dtype=object))
            blast.loc[red_index, 'homo' + str(i)] = 1
            blast.loc[blue_index, 'homo' + str(i)] = 0
            blast.loc[gray_index, 'homo' + str(i)] = -1
        blast.index = blast[0] + ',' + blast[1]
        return blast

    # AnchorWave pro collinearity file features :left gene is ref and right gene is query and ref_chr&query_chr
    def run(self):
        block_list = []   # collinearity_gene_pair to calculate tandem_ratio
        index_list = []   # block number, 1,2,3,4,5, ... ,num(block)
        chr1 = []         # query chr
        chr2 = []         # ref chr

        start1 = []       # query first line, left and top       start1 > end1 or start1 < end1
        end1 = []         # query last line, left and bottom
        start2 = []       # ref first line, right and top        start2 > end2
        end2 = []         # ref last line, right and bottom
        block_len = []    # element is block_length, length equal to num(block)

        block_order = []  # all_order object's element, every block's two seq order

        block_ks = []     # all_order object's element, every block's element ks
        ks_median_list = []
        ks_average_list = []

        block1 = []       # element is str which is block query order concatenated by "_"
        block2 = []       # element is str which is block ref order concatenated by "_"
        ks = []           # element is str which is block's every ks_value concatenated by "_"
        tandem_list = []
        block_homo = []
        homo_list = []
        ks_series = read_ks(self.ks, self.ks_col)
        table_df, map2 = read_table_file(self.table_file)
        blast = self.blast_homo(self.blast_file, self.repeat_number)
        # print("line139")
        with (open(self.collinearity) as f):
            _ = next(f)
            _ = next(f)
            for line in f:
                if line.startswith("#"):
                    if block_order and block_ks:
                        # process second block and record the first block's part info
                        start1, start2, end1, end2, block_order = record_block_region(start1, start2, end1, end2, block_order)
                        ks_median = get_median(block_ks)
                        ks_average = get_average(block_ks)
                        ks_median_list.append(ks_median)
                        ks_average_list.append(ks_average)
                        block1, block2, ks, block_ks, block_order \
                            = use_underscore_concatenate(block1, block2, ks, block_ks, block_order)
                        block = pd.DataFrame(block_list, columns=['id1', 'id2'])
                        tandem_list.append(self.tandem_ratio(table_df, map2, block))
                        df = pd.DataFrame(block_homo)
                        homo = df.mean().values.tolist()
                        if len(homo) == 0:
                            homo = [-1, -1, -1, -1, -1]
                        homo_list.append(homo)
                        block_list = []
                        block_order = []
                        block_ks = []
                        block_homo = []
                    # process the first section and the first block
                    index_list.append(int(line.split()[1]))
                    block_len.append(int(line.split()[2].split(sep="N=")[1]))
                    chr_pair = line.split()[4].split(sep="&")
                    chr1.append(str(chr_pair[1]))
                    chr2.append(str(chr_pair[0]))
                else:
                    # process the second section and the first block
                    split_line = line.split()
                    block_list.append((split_line[5], split_line[0]))
                    if str(split_line[5]) + "," + str(split_line[0]) in blast.index:
                        block_homo.append(blast.loc[str(split_line[5]) + ","
                                                    + str(split_line[0]), ['homo' + str(i) for i in range(1, 6)]].values.tolist())
                    if split_line[5] + "," + split_line[0] in ks_series.index:
                        ks_value = ks_series[split_line[5] + "," + split_line[0]]
                    elif split_line[0] + "," + split_line[5] in ks_series.index:
                        ks_value = ks_series[split_line[0] + "," + split_line[5]]
                    else:
                        ks_value = -1
                    block_ks.append(ks_value)
                    block_order.append((int(split_line[7]), int(split_line[2])))
            # the last block's part info
            start1, start2, end1, end2, block_order = record_block_region(start1, start2, end1, end2, block_order)
            ks_median = get_median(block_ks)
            ks_average = get_average(block_ks)
            ks_median_list.append(ks_median)
            ks_average_list.append(ks_average)
            block1, block2, ks, block_ks, block_order \
                = use_underscore_concatenate(block1, block2, ks, block_ks, block_order)
            block = pd.DataFrame(block_list, columns=['id1', 'id2'])
            tandem_list.append(self.tandem_ratio(table_df, map2, block))
            df = pd.DataFrame(block_homo)
            homo = df.mean().values.tolist()
            if len(homo) == 0:
                homo = [-1, -1, -1, -1, -1]
            homo_list.append(homo)
            [homo1, homo2, homo3, homo4, homo5] = list(zip(*homo_list))
            homo1 = list(homo1)
            homo2 = list(homo2)
            homo3 = list(homo3)
            homo4 = list(homo4)
            homo5 = list(homo5)

        data = [index_list, chr1, chr2, start1, end1, start2, end2, block_len, ks_median, ks_average, homo1, homo2,
                homo3, homo4, homo5, block1, block2, ks, tandem_list]
        columns = ['id', 'chr1', 'chr2', 'start1', 'end1', 'start2', 'end2', 'length', 'ks_median', 'ks_average',
                   'homo1', 'homo2', 'homo3', 'homo4', 'homo5', 'block1', 'block2', 'ks', 'tandem_ratio']
        data_dict = dict(zip(columns, data))
        df = pd.DataFrame(data_dict)
        df['density1'] = df['length'] / ((df['end1'] - df['start1']).abs() + 1)
        df['density2'] = df['length'] / ((df['end2'] - df['start2']).abs() + 1)
        df['class1'] = 0
        df['class2'] = 0
        df.to_csv(self.savefile, index=False)
