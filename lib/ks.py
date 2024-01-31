class Ks:
    def __init__(self):
        self.refpepfile = ''
        self.qurpepfile = ''
        self.refcdsfile = ''
        self.qurcdsfile = ''
        self.collinearityfile = ''

    def run_ks(self):
        pairs = []
        refpep_dict = {}
        refcds_dict = {}
        qurpep_dict = {}
        qurcds_dict = {}
        with open('coll', 'r') as file:
            for line_number, line in enumerate(file, 1):
                line = line.strip()
                if line_number == 2 or line.startswith('#'):
                    continue
                fields = line.split()
                if len(fields) >= 6:
                    gene_name_1 = fields[0]
                    gene_name_2 = fields[5]
                    pairs.append(gene_name_1)
                    pairs.append(gene_name_2)

                with open('refpep', 'r') as pep_file:
                    current_gene_name = None
                    current_gene_sequence = ""
                    for line in pep_file:
                        line = line.strip()
                        if line.startswith('>'):
                            current_gene_name = line[1:]
                        else:
                            if current_gene_name in pairs:
                                current_gene_sequence += line
                                if current_gene_name and current_gene_name in pairs:
                                    refpep_dict[current_gene_name] = current_gene_sequence
                                    print(refpep_dict)
                                    with open('refpepfile', 'w') as file:
                                        for key, value in refpep_dict.items():
                                            file.write(f">{key}\n{value}\n")

                with open('qurpep', 'r') as pep_file:
                    current_gene_name = None
                    current_gene_sequence = ""
                    for line in pep_file:
                        line = line.strip()
                        if line.startswith('>'):
                            current_gene_name = line[1:]
                        else:
                            if current_gene_name in pairs:
                                current_gene_sequence += line
                                if current_gene_name and current_gene_name in pairs:
                                    qurpep_dict[current_gene_name] = current_gene_sequence
                                    print(qurpep_dict)
                                    with open('qurpepfile', 'w') as file:
                                        for key, value in qurpep_dict.items():
                                            file.write(f">{key}\n{value}\n")

                with open('refcds', 'r') as cds_file:
                    current_gene_name = None
                    current_gene_sequence = ""
                    for line in cds_file:
                        line = line.strip()
                        if line.startswith('>'):
                            current_gene_name = line[1:]
                        else:
                            if current_gene_name in pairs:
                                current_gene_sequence += line
                                if current_gene_name and current_gene_name in pairs:
                                    refcds_dict[current_gene_name] = current_gene_sequence
                                    print(refcds_dict)
                                    with open('refcdsfile', 'w') as file:
                                        for key, value in refcds_dict.items():
                                            file.write(f">{key}\n{value}\n")

                with open('qurpep', 'r') as pep_file:
                    current_gene_name = None
                    current_gene_sequence = ""
                    for line in pep_file:
                        line = line.strip()
                        if line.startswith('>'):
                            current_gene_name = line[1:]
                        else:
                            if current_gene_name in pairs:
                                current_gene_sequence += line
                                if current_gene_name and current_gene_name in pairs:
                                    qurpep_dict[current_gene_name] = current_gene_sequence
                                    print(qurpep_dict)
                                    with open('qurpepfile', 'w') as file:
                                        for key, value in qurpep_dict.items():
                                            file.write(f">{key}\n{value}\n")

                with open('qurcds', 'r') as cds_file:
                    current_gene_name = None
                    current_gene_sequence = ""
                    for line in cds_file:
                        line = line.strip()
                        if line.startswith('>'):
                            current_gene_name = line[1:]
                        else:
                            if current_gene_name in pairs:
                                current_gene_sequence += line
                                if current_gene_name and current_gene_name in pairs:
                                    qurcds_dict[current_gene_name] = current_gene_sequence
                                    print(qurcds_dict)
                                    with open('qurcdsfile', 'w') as file:
                                        for key, value in qurcds_dict.items():
                                            file.write(f">{key}\n{value}\n")

                refcds_dict = {}
                refpep_dict = {}
                qurcds_dict = {}
                qurpep_dict = {}
                pairs = []


test = Ks()
test.run_ks()