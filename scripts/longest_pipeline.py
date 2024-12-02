import subprocess
import os, sys, argparse, logging, re
from argparse import ArgumentParser


logger = logging.getLogger('longest')
logger.setLevel(level=logging.INFO)
logging.basicConfig(level=logging.INFO,
                    datefmt='%Y/%m/%d %H:%M:%S',
                    format='[%(asctime)s %(levelname)s] %(message)s')

SOFTWARE = "quota_Anchor"
RAW_PEP_SUFFIX = ".raw.pep"
LONGEST_PEP_SUFFIX = ".longest.pep"
RAW_CDS_SUFFIX = ".raw.cds"
LONGEST_CDS_SUFFIX = ".longest.cds"
GFF3_SUFFIX = ".gff3"
FASTA_SUFFIX = ".fa"


def search_input_dir_files(input_directory):
    _genome_lst = []
    _gff_lst = []
    _species_lst = []

    abs_input_dir = os.path.abspath(input_directory)

    files_and_folders = os.listdir(abs_input_dir)
    for name in files_and_folders:
        abs_path = os.path.join(abs_input_dir, name)
        if os.path.isfile(abs_path) and name.endswith(FASTA_SUFFIX):
            _genome_lst.append(abs_path)
            _species_lst.append(name[:-3])
        elif os.path.isfile(abs_path) and name.endswith(GFF3_SUFFIX):
            _gff_lst.append(abs_path)
        else:
            pass
    _genome_lst.sort()
    _gff_lst.sort()
    _species_lst.sort()
    return _genome_lst, _gff_lst, _species_lst

def output_info(info):
    if info:
        info_list = info.split("\n")
        for i in info_list:
            if i:
               print(f'{i}')
    print()
    return

class Longest:
    def __init__(self, overwrite):
        self.overwrite = overwrite

    @staticmethod
    def get_longest_four_parameters(genome_lst, gff_lst, species_lst, output_path):
        genome_str = ",".join(genome_lst)
        gff_str = ",".join(gff_lst)

        raw_pep_str = ""
        longest_pep_str = ""
        raw_cds_str = ""
        longest_cds_str = ""

        abs_output_path = os.path.abspath(output_path)
        longest_path = os.path.join(abs_output_path, "01longest")
        os.makedirs(longest_path, exist_ok=True)

        for i in species_lst:
            raw_pep_str =  raw_pep_str + os.path.join(longest_path, str(i) + RAW_PEP_SUFFIX) + ","
            longest_pep_str = longest_pep_str + os.path.join(longest_path, str(i) + LONGEST_PEP_SUFFIX) + ","
            raw_cds_str = raw_cds_str + os.path.join(longest_path, str(i) + RAW_CDS_SUFFIX) + ","
            longest_cds_str = longest_cds_str + os.path.join(longest_path, str(i) + LONGEST_CDS_SUFFIX) + ","
        merged_pep = os.path.join(longest_path, "merged.longest.pep")
        merged_cds = os.path.join(longest_path, "merged.longest.cds")
        return genome_str, gff_str, raw_pep_str, longest_pep_str, raw_cds_str, longest_cds_str, merged_pep, merged_cds

    def quota_anchor_longest_pep(self, fasta, gff3, raw_pep_var, longest_pep_var, merged_pep_var, species_number):
        if self.overwrite:
            command_line = [SOFTWARE,
                            'longest_pep',
                            '-f', fasta,
                            '-g', gff3,
                            '-p', raw_pep_var,
                            '-l', longest_pep_var,
                            '-t', species_number,
                            '-merge', merged_pep_var,
                            '--overwrite'
                            ]
        else:
            command_line = [SOFTWARE,
                            'longest_pep',
                            '-f', fasta,
                            '-g', gff3,
                            '-p', raw_pep_var,
                            '-l', longest_pep_var,
                            '-merge', merged_pep_var,
                            '-t', species_number
                            ]
        try:
            # For debugging, output will be delayed
            result = subprocess.run(command_line, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
            stderr_pep = result.stderr
            stdout_pep = result.stdout
            output_info(stderr_pep)
            output_info(stdout_pep)

        except subprocess.CalledProcessError as e:
            pattern = re.compile(r'(set)\s+\'(--overwrite)\'\s+(in)\s+(the)\s+(command)\s+(line)\s+(to)\s+(overwrite)\s+(it)\.')
            error_message = e.stderr
            output_info(error_message)

            output_message = e.stdout
            output_info(output_message)

            if not re.search(pattern, error_message) and not re.search(pattern, output_message):
                sys.exit(1)
            pass

    def quota_anchor_longest_cds(self, fasta, gff3, raw_cds_var, longest_cds_var, merged_cds_var, species_number):
        if self.overwrite:
            command_line = [SOFTWARE,
                            'longest_cds',
                            '-f', fasta,
                            '-g', gff3,
                            '-p', raw_cds_var,
                            '-l', longest_cds_var,
                            '-merge', merged_cds_var,
                            '-t', species_number,
                            '--overwrite'
                            ]
        else:
            command_line = [SOFTWARE,
                            'longest_cds',
                            '-f', fasta,
                            '-g', gff3,
                            '-p', raw_cds_var,
                            '-l', longest_cds_var,
                            '-merge', merged_cds_var,
                            '-t', species_number
                            ]
        try:
            result = subprocess.run(command_line, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
            stderr_cds = result.stderr
            stdout_cds = result.stdout
            output_info(stderr_cds)
            output_info(stdout_cds)

        except subprocess.CalledProcessError as e:
            pattern = re.compile(r'(set)\s+\'(--overwrite)\'\s+(in)\s+(the)\s+(command)\s+(line)\s+(to)\s+(overwrite)\s+(it)\.')
            error_message = e.stderr
            output_info(error_message)

            output_message = e.stdout
            output_info(output_message)

            if not re.search(pattern, error_message) and not re.search(pattern, output_message):
                sys.exit(1)
            pass

    def run(self):
        genome_files, gff_files, species_lst = search_input_dir_files(INPUT_DIR)
        parallel_number = str(len(species_lst))
        genome_str, gff_str, raw_pep_str, longest_pep_str, raw_cds_str, longest_cds_str, merged_pep, merged_cds = self.get_longest_four_parameters(genome_files, gff_files, species_lst, OUTPUT_DIR)
        if not SKIP_LONGEST_PEP:
            self.quota_anchor_longest_pep(genome_str, gff_str, raw_pep_str, longest_pep_str, merged_pep, parallel_number)
        if not SKIP_LONGEST_CDS:
            self.quota_anchor_longest_cds(genome_str, gff_str, raw_cds_str, longest_cds_str, merged_cds, parallel_number)
        return species_lst


if __name__ == '__main__':
    parser = ArgumentParser(usage="Generate longest protein and longest cds for each species in the input directory.", formatter_class=argparse.RawDescriptionHelpFormatter, description=
    """
    The initial goal of this script is to simplify the user pipeline of quota_Anchor positioning wgd events relative to species divergent events.
    Note:
    1. The directory corresponding to the input_dir parameter needs to contain genome files and genome annotation files with .gff3 suffix and .fa suffix.
    
    2. The following is the current directory tree.
                        ├── raw_data
                        │   ├── maize.fa
                        │   ├── maize.gff3
                        │   ├── oryza.fa
                        │   ├── oryza.gff3
                        │   ├── setaria.fa
                        │   ├── setaria.gff3
                        │   ├── sorghum.fa
                        │   └── sorghum.gff3
                        └── scripts
                            ├── ks_pipeline.py
                            └── longest_pipeline.py
    
    3. The script will create input_dir/01longest folders in current directory.
                    ├── output_dir
                    │     └── 01longest
                    ├── raw_data
                    │     ├── maize.fa
                    │     ├── maize.gff3
                    │     ├── oryza.fa
                    │     ├── oryza.gff3
                    │     ├── setaria.fa
                    │     ├── setaria.gff3
                    │     ├── sorghum.fa
                    │     └── sorghum.gff3
                    ├── scripts
                        ├── ks_pipeline.py
                        └── longest_pipeline.py
    4. In other words, You need to provide genome files in fasta format and annotation files in gff3 format.
    
    Example:
    1. python ./scripts/longest_pipeline.py -i raw_data -o output_dir
    2. python ./scripts/longest_pipeline.py -i raw_data -o output_dir --overwrite
    2. python ./scripts/longest_pipeline.py -i raw_data -o output_dir --overwrite --skip_longest_pep
    3. python ./scripts/longest_pipeline.py -i raw_data -o output_dir --overwrite --skip_longest_cds
    """)
    parser.add_argument("-i", "--input_dir",
                        dest="input_dir",
                        type=str,
                        default="",
                        metavar="",
                        help="Path to the directory containing genome files and genome annotation files with .gff3 and .fa suffix")
    parser.add_argument("-o", "--output_dir",
                        dest="output_dir",
                        type=str,
                        default="",
                        metavar="",
                        help="If the path does not exist, the script will recursively create folders in this path.")
    parser.add_argument("-skip_longest_pep", "--skip_longest_pep",
                        dest="skip_longest_pep",
                        action='store_true',
                        help="Skip generating longest cds files for each species in the input directory.")
    parser.add_argument("-skip_longest_cds", "--skip_longest_cds",
                        dest="skip_longest_cds",
                        action='store_true',
                        help="Skip generating longest cds files for each species in the input directory.")
    parser.add_argument("-overwrite", "--overwrite",
                        dest="overwrite",
                        action='store_true',
                        help="Overwrite the output file.")
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    if args.input_dir == "":
        logger.error("Error: please specify --input_dir.")
        parser.print_help()
        sys.exit(1)
    if args.output_dir == "":
        logger.error("Error: please specify --output_dir.")
        parser.print_help()
        sys.exit(1)

    INPUT_DIR = args.input_dir
    OUTPUT_DIR = args.output_dir
    SKIP_LONGEST_PEP = args.skip_longest_pep
    SKIP_LONGEST_CDS = args.skip_longest_cds
    OVERWRITE = args.overwrite

    SPECIES_LST = Longest(OVERWRITE).run()
