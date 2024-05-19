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

    def run_coll_dotplot(self):
        dict2 = {"color": "strand"}
        query_chr = self.special_query_chr.split(sep=",")
        ref_chr = self.special_ref_chr.split(sep=",")
        coll_df = pd.read_csv(self.input_file, header=0, sep="\t", comment="#", low_memory=False)
        coll_df["queryChr"] = coll_df["queryChr"].astype(str)
        coll_df["refChr"] = coll_df["refChr"].astype(str)
        coll_df = coll_df[coll_df['queryChr'].isin(query_chr)]
        coll_df = coll_df[coll_df['refChr'].isin(ref_chr)]
        if self.type == 'order':
            coll_df['queryId'] = coll_df['queryId'].apply(lambda x: x / 1000)
            coll_df['refId'] = coll_df['refId'].apply(lambda x: x / 1000)
            dict1 = {"x": "queryId", "y": "refId"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=1) +
                    scale_x_continuous() +
                    scale_y_continuous())

            my_theme = theme(
                axis_line=element_blank(),
                panel_background=element_rect(fill='white'),
                panel_border=element_rect(fill=None, color="black", linewidth=2, linetype="solid"),
                # axis_text_y=element_text(color="black"),
                legend_position=element_blank(),
                axis_ticks_major_x=True,
                axis_ticks_major_y=True,
                axis_title_x=str(self.query_name) + "(unit:k)",
                axis_title_y=str(self.ref_name) + "(unit:k)",
                axis_text=element_text(size=2),
                axis_text_x=element_text(rotation=200.0, ha='center', va='bottom', color="red")
            )

            plot1 = plot + my_theme + theme_grey(base_size=40) + labs(x=f'{self.query_name} (Unit:K)', y=f'{self.ref_name} (Unit:K)')
            plot1.save(str(self.filename), width=1500, height=1200, units="mm", limitsize=False)
        else:
            coll_df['queryStart'] = coll_df['queryStart'].apply(lambda x: x / 1000000)
            coll_df['referenceStart'] = coll_df['referenceStart'].apply(lambda x: x / 1000000)
            dict1 = {"x": "queryStart", "y": "referenceStart"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=1) +
                    scale_x_continuous() +
                    scale_y_continuous())

            my_theme = theme(
                axis_line=element_blank(),
                panel_background=element_rect(fill='white'),
                panel_border=element_rect(fill=None, color="black", linewidth=2, linetype="solid"),
                axis_text_y=element_text(color="black"),
                legend_position=element_blank(),
                axis_ticks_major_x=True,
                axis_ticks_major_y=True,
                axis_title_x=str(self.query_name) + "(unit:k)",
                axis_title_y=str(self.ref_name) + "(unit:k)",
                axis_text=element_text(size=2),
                axis_text_x=element_text(rotation=200.0, ha='center', va='bottom', color="red")
            )

            plot1 = plot + my_theme + theme_grey(base_size=40) + labs(x=f'{self.query_name} (Unit:M)', y=f'{self.ref_name} (Unit:M)')
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
        if self.type == "order":
            dict1 = {"x": "queryId", "y": "refId"}
            coll_df['refId'] = coll_df['refId'].apply(lambda x: x / 1000)
            coll_df['queryId'] = coll_df['queryId'].apply(lambda x: x / 1000)

            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=1) +
                    scale_x_continuous() +
                    scale_y_continuous())

            my_theme = theme(
                axis_line=element_blank(),
                panel_background=element_rect(fill='white'),
                panel_border=element_rect(fill=None, color="black", linewidth=2, linetype="solid"),
                axis_text_y=element_text(color="black"),
                legend_position=element_blank(),
                axis_ticks_major_x=True,
                axis_ticks_major_y=True,
                # axis_title_x=str(self.query_name) + "(unit:K)",
                # axis_title_y=str(self.ref_name) + "(unit:K)",
                axis_text=element_text(size=2),
                axis_text_x=element_text(rotation=200.0, ha='center', va='bottom', color="red")
            )

            plot1 = plot + my_theme + theme_grey(base_size=40) + labs(x=f'{self.query_name} (Unit:K)', y=f'{self.ref_name} (Unit:M)')
            plot1.save(str(self.filename), width=1500, height=1200, units="mm", limitsize=False)
        else:
            coll_df['queryStart'] = coll_df['queryStart'].apply(lambda x: x / 1000000)
            coll_df['referenceStart'] = coll_df['referenceStart'].apply(lambda x: x / 1000000)
            dict1 = {"x": "queryStart", "y": "referenceStart"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=1) +
                    scale_x_continuous() +
                    scale_y_continuous())

            my_theme = theme(
                axis_line=element_blank(),
                panel_background=element_rect(fill='white'),
                panel_border=element_rect(fill=None, color="black", linewidth=2, linetype="solid"),
                axis_text_y=element_text(color="black"),
                legend_position=element_blank(),
                axis_ticks_major_x=True,
                axis_ticks_major_y=True,
                axis_title_x=str(self.query_name) + "(unit:M)",
                axis_title_y=str(self.ref_name) + "(unit:M)",
                axis_text=element_text(size=2),
                axis_text_x=element_text(rotation=200.0, ha='center', va='bottom', color="red")
            )

            plot1 = plot + my_theme + theme_grey(base_size=40) + labs(x=f'{self.query_name} (Unit:M)', y=f'{self.ref_name} (Unit:M)')
            plot1.save(str(self.filename), width=1500, height=1200, units="mm", limitsize=False)

    def run(self):
        with open(self.input_file, 'r') as file:
            first_line = file.readline()
        if first_line.startswith("#"):
            self.run_coll_dotplot()
        else:
            self.run_blast_dotplot()
