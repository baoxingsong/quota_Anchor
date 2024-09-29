import os
import sys


class FileEmptyError(Exception):
    pass


def split_conf(conf, separator):
    new_conf = []
    split_lt = conf.split(separator)
    for ele in split_lt:
        ele = ele.strip()
        if len(ele) == 0:
            continue
        new_conf.append(ele)
    return new_conf


def read_collinearity(qry_prefix, ref_prefix, collinearity, chr_list, chr_to_start):
    # left -> ref;  right -> query
    data = []
    ref_chr_list = []
    query_chr_list = []
    gene_pos_dict = {}
    block_index = 0
    block = []
    # print(qry_prefix)
    # print(ref_prefix)
    # print(collinearity)
    # print(chr_list)
    # print(chr_to_start)
    with open(collinearity) as f:
        print("read", collinearity, "....")
        _ = next(f)
        _ = next(f)
        flag = True
        for line in f:
            if line.startswith("#"):
                if block:
                    data.append([block[0][0], block[0][1], block[-1][0], block[-1][1]])
                    block = []
                chr_pair = line.split()[4].split("&")
                if ref_prefix + chr_pair[0] not in chr_list or qry_prefix + chr_pair[1] not in chr_list:
                    flag = False
                else:
                    flag = True
                    ref_chr = ref_prefix + chr_pair[0]
                    ref_chr_list.append(ref_chr)
                    query_chr = qry_prefix + chr_pair[1]
                    query_chr_list.append(query_chr)
                    block_index += 1
            else:
                if flag:
                    line_list = line.split()
                    block.append([line_list[0], line_list[5]])
                    gene_pos_dict[line_list[0]] = chr_to_start[ref_chr] + int(line_list[4])
                    gene_pos_dict[line_list[5]] = chr_to_start[query_chr] + int(line_list[9])
                else:
                    continue
        data.append([block[0][0], block[0][1], block[-1][0], block[-1][1]])
        print("parse", collinearity, "success")
    return data, gene_pos_dict, ref_chr_list, query_chr_list


def bezier3(a, b, c, d, t):
    fvalue = []
    for i in t:
        y = a * (1-i)**3 + 3 * b * i * (1-i)**2 + 3 * c * (1-i) * i**2 + d * i**3
        fvalue.append(y)
    return fvalue


def file_empty(file_path):
    try:
        file_path = os.path.abspath(file_path)
        if os.path.exists(file_path):
            if os.path.getsize(file_path) > 0:
                pass
            else:
                raise FileEmptyError
        else:
            raise FileNotFoundError
    except FileNotFoundError as e1:
        exist_error_message = f"{file_path} don't exist or the path physically exists but permission is not granted to execute os.stat() on the requested file"
        print(exist_error_message)  
        sys.exit(1)
    except FileEmptyError as e2:
        empty_error_message = "{0} is empty".format(file_path)
        print(empty_error_message)
        sys.exit(1)
    except OSError as e3:
        OSE_error_message = f"{file_path} does not exist or is inaccessible."
        print(OSE_error_message)
        sys.exit(1)


def output_file_parentdir_exist(path, overwite):
    if os.path.exists(path):
        print(f"{path} already exist.")
        if overwite:
            print(f"{path} will be overwrited.")
            return
        else:
            print(f"{path} will not be overwrited, and you can set '--overwrite' in the command line to overwrite it.")
            sys.exit(1)
    path = os.path.abspath(path)
    dir_name = os.path.dirname(path)
    if os.path.isdir(dir_name):
        pass
    else:
        print(f"{dir_name} does not exist and software will make directorys.")
        os.makedirs(dir_name, exist_ok=True)

