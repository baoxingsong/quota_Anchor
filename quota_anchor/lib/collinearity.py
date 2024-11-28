import subprocess
import resource
import sys
from . import base
import logging

logger = logging.getLogger('main.collinearity')

class Collinearity:
    def __init__(self, config_pra, config_soft, parameter):
        self.overwrite = False
        self.r_value = 1
        self.q_value = 1
        self.maximum_gap_size = 25
        self.tandem_length = 0
        self.over_lap_window = 1
        self.count_style = 0
        self.strict_strand = 1
        self.get_all_collinearity = 0
        self.minimum_chain_score = 3
        self.gap_extend_penalty = -0.005
        self.strict_remove_overlap = 0
        self.AnchorWave = config_soft['software']['AnchorWave']

        self.input_file_name = ""
        self.output_coll_name = ""
        self.input_file_name = ""
        self.output_coll_name =""

        for i in config_pra.sections():
            if i == 'AnchorWave':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)

    def convert_attr_type(self):
        attr = ['r_value', 'q_value', 'maximum_gap_size', 'tandem_length', 'over_lap_window', 'count_style',
                'strict_strand', 'get_all_collinearity',  'minimum_chain_score', 'gap_extend_penalty', 'strict_remove_overlap']
        for i in attr:
            setattr(self, i, str(getattr(self, i)))

    def run_anchorwave_pro(self,
                        input_file,
                        output_file,
                        r_value,
                        q_value,
                        maximum_gap_size,
                        tandem_length,
                        over_lap_window,
                        count_style,
                        strict_strand,
                        get_all_collinearity,
                        minimum_chain_score,
                        gap_extend_penalty,
                        strict_remove_overlap
                        ):
        command_line = [self.AnchorWave,
                        'pro', 
                        '-i', input_file,
                        '-o', output_file,
                        '-R', r_value,
                        '-Q', q_value,
                        '-D', maximum_gap_size,
                        '-m', tandem_length,
                        '-W', over_lap_window,
                        '-c', count_style,
                        '-s', strict_strand,
                        '-a', get_all_collinearity,
                        '-I', minimum_chain_score,
                        '-E', gap_extend_penalty,
                        '-f', strict_remove_overlap
                        ]
        try:
            logger.info(f"Run AnchorWave and generate {output_file} start.")
            result = subprocess.run(command_line, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_gff_read = result.stderr.decode()
            stdout_gff_read = result.stdout.decode()
            base.output_info(stderr_gff_read)
            base.output_info(stdout_gff_read)
            logger.info(f"Run AnchorWave and generate {output_file} end.")
        # Added in version 3.5
        except subprocess.CalledProcessError as e:
            logger.error(f"Run AnchorWave and generate {output_file} failed!")

            error_message = e.stderr.decode()
            base.output_info(error_message)

            output_message = e.stdout.decode()
            base.output_info(output_message)
            sys.exit(1)

    def col_init(self):
        if not self.input_file_name:
            logger.error("Please specify your table file as input file")
            sys.exit(1)

        if not self.output_coll_name:
            logger.error("Please specify your output file name")
            sys.exit(1)


    def run(self):
        logger.info("Collinearity module init and the following parameters are config information.")
        print()
        for key, value in vars(self).items():
            if key != "AnchorWave" and key != "conf":
                print(key, "=", value)
        print()
        self.col_init()
        base.file_empty(self.input_file_name)
        base.output_file_parentdir_exist(self.output_coll_name, self.overwrite)

        # stack_size = 1024 * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, -1))

        self.convert_attr_type()
        self.run_anchorwave_pro(self.input_file_name,
                                self.output_coll_name, 
                                self.r_value,
                                self.q_value,
                                self.maximum_gap_size,
                                self.tandem_length,
                                self.over_lap_window,
                                self.count_style,
                                self.strict_strand,
                                self.get_all_collinearity,
                                self.minimum_chain_score,
                                self.gap_extend_penalty,
                                self.strict_remove_overlap
                                )
        logger.info(f"Collinearity analysis process finished!")
