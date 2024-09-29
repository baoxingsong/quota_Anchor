import subprocess
from . import base
import os 


class Collinearity:
    def __init__(self, config_pra, config_soft, parameter):
        
        self.overwrite = False
        self.r_value = str(1)
        self.q_value = str(1)
        self.maximum_gap_size = str(25)
        self.tandem_length = str(0)
        self.over_lap_window = str(0)
        self.count_style = str(1)
        self.strict_strand = str(1)
        self.get_all_collinearity = str(0)
        self.minimum_chain_score = str(2)
        self.AnchorWave = config_soft['software']['AnchorWave']
        
        for i in config_pra.sections():
            if i == 'AnchorWave':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)
        print()
        for key, value in vars(self).items():
            if key != "AnchorWave" and key != "conf":
                print(key, "=", value)
        print()


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
                        minimum_chain_score
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
                        '-I', minimum_chain_score
                        ]
        try:
            result = subprocess.run(command_line, check=True)
        # Added in version 3.5
        except subprocess.CalledProcessError:
            print(f"{result}, class Collinearity's function run_anchorwave_pro failed, and return a non-zero code.")

    def run(self):
        base.file_empty(self.input_file_name)
        base.output_file_parentdir_exist(self.output_coll_name, self.overwrite)
        self.run_anchorwave_pro(self.input_file_name,
                                self.output_coll_name, 
                                str(self.r_value), 
                                str(self.q_value),
                                str(self.maximum_gap_size), 
                                str(self.tandem_length), 
                                str(self.over_lap_window), 
                                str(self.count_style), 
                                str(self.strict_strand), 
                                str(self.get_all_collinearity),
                                str(self.minimum_chain_score))
        return self.output_coll_name
