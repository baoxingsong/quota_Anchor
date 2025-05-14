import sys


ks_par = ""
col_par = ""
with open(sys.argv[1], 'r') as f:
    _ = next(f)
    for line in f:
        line_list = line.split(",")
        ks_file =  "./output_dir/02synteny/" + line_list[0] + "_" +  line_list[1] + "0" + ".ks,"
        collinearity_file =  "./output_dir/02synteny/" + line_list[0] + "_" +  line_list[1] + "0" + ".collinearity,"
        ks_par += ks_file
        col_par += collinearity_file


print(ks_par)
print(col_par)