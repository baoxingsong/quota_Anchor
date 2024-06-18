from plotnine import *
import pandas as pd
import numpy as np


class Dotplot:
    def __init__(self, config_pra):
        self.input_file = config_pra['dotplot']['input_file']
        self.special_query_chr = config_pra['dotplot']['special_query_chr']
        self.special_ref_chr = config_pra['dotplot']['special_ref_chr']
        self.type = config_pra['dotplot']['type']
        self.query_name = config_pra['dotplot']['query_name']
        self.ref_name = config_pra['dotplot']['ref_name']
        self.filename = config_pra['dotplot']['filename']
        self.my_theme = theme(
                panel_background=element_rect(fill='white'),
                panel_border=element_rect(fill=None, color="black", linewidth=2, linetype="solid"),
                axis_ticks_major_x=element_line(size=4),
                axis_ticks_major_y=element_line(size=4),
                legend_text=element_text(size=50, ha='center'),
                legend_title=element_text(size=50, ha='center'),
                axis_text_x=element_text(rotation=-60, hjust=0, vjust=1)
            )

    @staticmethod
    def major_formatter(breaks):
        return [f'{int(b)}M' for b in breaks]

    def run_coll_dotplot(self):
        dict2 = {"color": "strand"}
        query_chr = self.special_query_chr.split(sep=",")
        ref_chr = self.special_ref_chr.split(sep=",")
        coll_df = pd.read_csv(self.input_file, header=0, sep="\t", comment="#", low_memory=False)
        coll_df["queryChr"] = coll_df["queryChr"].astype(str)
        coll_df["refChr"] = coll_df["refChr"].astype(str)
        coll_df = coll_df[coll_df['queryChr'].isin(query_chr)]
        coll_df = coll_df[coll_df['refChr'].isin(ref_chr)]
        custom_colors = {"+": "red", "-": "blue"}
        if self.type == 'order':
            dict1 = {"x": "queryId", "y": "refId"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=2, alpha=1) +
                    scale_color_manual(values=custom_colors) +
                    scale_x_continuous() +
                    scale_y_continuous())

            plot1 = plot + theme_grey(base_size=50) + labs(x=f'{self.query_name}', y=f'{self.ref_name}')
            plot1 = plot1 + guides(**{'color': guide_legend(override_aes={'size': 8, 'alpha': 1})}) + self.my_theme
            plot1.save(str(self.filename), width=1500, height=1200, units="mm", limitsize=False)
        else:
            coll_df['queryStart'] = coll_df['queryStart'].apply(lambda x: x / 1000000)
            coll_df['referenceStart'] = coll_df['referenceStart'].apply(lambda x: x / 1000000)
            dict1 = {"x": "queryStart", "y": "referenceStart"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=2, alpha=1) +
                    scale_color_manual(values=custom_colors) +
                    scale_x_continuous(labels=self.major_formatter) +
                    scale_y_continuous(labels=self.major_formatter))

            plot1 = plot + theme_grey(base_size=50) + self.my_theme + labs(x=f'{self.query_name}', y=f'{self.ref_name}')
            plot1 = plot1 + guides(**{'color': guide_legend(override_aes={'size': 8, 'alpha': 1})})
            plot1.save(str(self.filename), width=1500, height=1200, units="mm", limitsize=False)

    def run_blast_dotplot(self):
        dict2 = {"color": "strand"}
        query_chr = self.special_query_chr.split(sep=",")
        ref_chr = self.special_ref_chr.split(sep=",")
        coll_df = pd.read_csv(self.input_file, header=None, sep="\t", comment="#", low_memory=False)
        coll_df.columns = ["refGene", "refChr", "refId", "referenceStart", "referenceEnd", "refStrand",
                           "queryGene",	"queryChr", "queryId", "queryStart", "queryEnd", "queryStrand", "identity"]
        coll_df["strand"] = np.where(coll_df["refStrand"] == coll_df["queryStrand"], '+', '-')
        coll_df["queryChr"] = coll_df["queryChr"].astype(str)
        coll_df["refChr"] = coll_df["refChr"].astype(str)
        coll_df = coll_df[coll_df['queryChr'].isin(query_chr)]
        coll_df = coll_df[coll_df['refChr'].isin(ref_chr)]
        custom_colors = {"+": "red", "-": "blue"}
        if self.type == "order":
            dict1 = {"x": "queryId", "y": "refId"}
            # coll_df['refId'] = coll_df['refId'].apply(lambda x: x / 1000)
            # coll_df['queryId'] = coll_df['queryId'].apply(lambda x: x / 1000)

            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=1, alpha=0.5) +
                    scale_color_manual(values=custom_colors) +
                    scale_x_continuous() +
                    scale_y_continuous())

            plot1 = plot + theme_grey(base_size=50) + self.my_theme + labs(x=f'{self.query_name}', y=f'{self.ref_name}')
            plot1 = plot1 + guides(**{'color': guide_legend(override_aes={'size': 8, 'alpha': 1})})
            plot1.save(str(self.filename), width=1500, height=1200, units="mm", limitsize=False)
        else:
            coll_df['queryStart'] = coll_df['queryStart'].apply(lambda x: x / 1000000)
            coll_df['referenceStart'] = coll_df['referenceStart'].apply(lambda x: x / 1000000)
            dict1 = {"x": "queryStart", "y": "referenceStart"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=0.8, alpha=0.2) +
                    scale_color_manual(values=custom_colors) +
                    scale_x_continuous(labels=self.major_formatter) +
                    scale_y_continuous(labels=self.major_formatter))

            plot1 = plot + theme_grey(base_size=50) + self.my_theme + labs(x=f'{self.query_name}', y=f'{self.ref_name}')
            plot1 = plot1 + guides(**{'color': guide_legend(override_aes={'size': 8, 'alpha': 1})})
            plot1.save(str(self.filename), width=1500, height=1200, units="mm", limitsize=False)

    def run(self):
        with open(self.input_file, 'r') as file:
            first_line = file.readline()
        if first_line.startswith("#"):
            self.run_coll_dotplot()
        else:
            self.run_blast_dotplot()
