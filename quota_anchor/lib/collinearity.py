import subprocess


class Collinearity:
    def __init__(self, config_pra, config_soft):

        self.AnchorWave = config_soft['software']['AnchorWave']

        self.input_file_name = config_pra['AnchorWave']['input_file_name']
        self.R = config_pra["AnchorWave"]["R"]
        self.Q = config_pra["AnchorWave"]["Q"]
        self.maximum_gap_size = config_pra["AnchorWave"]["maximum_gap_size"]
        self.delete_tandem = config_pra["AnchorWave"]["delete_tandem"]
        self.tandem_dis = config_pra["AnchorWave"]["tandem_dis"]
        self.output_coll_name = config_pra["AnchorWave"]["output_coll_name"]

    def run_anchorwave_pro(self, input_file, output_file, r_value, q_value, maximum_gap_size, delete_tandem, tandem_dis):
        command_line = [self.AnchorWave, 'pro', '-i', input_file, '-o', output_file, '-R', r_value, '-Q', q_value, '-D', maximum_gap_size, '-m', delete_tandem, '-W', tandem_dis]
        try:
            subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError:
            print("class_Collinearity's function run_AnchorWave_pro failed:")

    def run(self):
        self.run_anchorwave_pro(self.input_file_name, self.output_coll_name, self.R, self.Q, self.maximum_gap_size, self.delete_tandem, self.tandem_dis)
        return self.output_coll_name
