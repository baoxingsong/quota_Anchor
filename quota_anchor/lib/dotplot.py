import logging
import sys

from plotnine import *
from plotnine.exceptions import PlotnineWarning
import pandas as pd
import numpy as np
from . import base
import matplotlib.pyplot as plt
from PIL import Image
import warnings
warnings.filterwarnings('ignore', category=PlotnineWarning)

logger = logging.getLogger('main.dotplot')
class Dotplot:
    def __init__(self, config_pra, parameter):
        self.overwrite = False
        self.italic = False
        self.type = "order"
        self.ref_name = "Reference species"
        self.query_name = "Query species"
        self.ks_area = "0,3"
        self.plotnine_figure_width=1500 
        self.plotnine_figure_height=1200
        self.disable_axis_text = False
        self.use_identity = False

        self.input_file = ""
        self.ref_length = ""
        self.query_length = ""
        self.output_file_name = ""
        self.ks = ""

        for i in config_pra.sections():
            if i == 'dotplot':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)

        self.col = "ks_NG86"
        if self.italic:
            self.my_theme = theme(
                    panel_background=element_rect(fill='white'),
                    panel_border=element_rect(fill=None, color="black", linewidth=2, linetype="solid"),
                    panel_spacing=0,

                    axis_title=element_text(size=90, face = "italic"),
                    
                    legend_text=element_text(size=40, ha='left', va='center'),
                    legend_title=element_text(size=40, ha='center', va='center'),
                    legend_ticks=element_line(size=4),

                    strip_text=element_text(size=50),
                    # strip_background=element_blank()
                    strip_background = element_rect(fill="white", color="black")
                ) 
        else:
            self.my_theme = theme(
                    panel_background=element_rect(fill='white'),
                    panel_border=element_rect(fill=None, color="black", linewidth=2, linetype="solid"),
                    panel_spacing=0,

                    axis_title=element_text(size=90),
                    
                    legend_text=element_text(size=40, ha='left', va='center'),
                    legend_title=element_text(size=40, ha='center', va='center'),
                    legend_ticks=element_line(size=4),

                    strip_text=element_text(size=50),
                    # strip_background=element_blank()
                    strip_background = element_rect(fill="white", color="black")
                )
    @staticmethod
    def major_formatter(breaks):
        return [f'{int(b)}M' for b in breaks]

    def prefix_rm(self, query_length, ref_length, df):
        query_df = self.read_length(query_length)
        ref_df = self.read_length(ref_length)
        if hasattr(self, "remove_chromosome_prefix"):
            prefix_list = self.remove_chromosome_prefix.split(",")
            for prf in prefix_list:
                df['queryChr'] = df['queryChr'].str.replace(rf'^{prf}', '', regex=True)
                df['refChr'] = df['refChr'].str.replace(rf'^{prf}', '', regex=True)
                query_df['chr'] = query_df['chr'].str.replace(rf'^{prf}', '', regex=True)
                ref_df['chr'] = ref_df['chr'].str.replace(rf'^{prf}', '', regex=True)
            return query_df, ref_df, df
        else:
            return query_df, ref_df, df

    @staticmethod
    def read_length(conf):
        df = pd.read_csv(conf, sep="\t", header=0, index_col=None)
        df['chr'] = df['chr'].astype(str)
        return df

    @staticmethod
    def read_ks(file, col):
        df = pd.read_csv(file, sep="\t", header=0, index_col=None, low_memory=False)
        df["id1"] = df["id1"].astype(str)
        df["id2"] = df["id2"].astype(str)

        df = df[["id1", "id2", col]]
        df_reverse = df[["id2", "id1", col]]
        df_reverse.rename(columns={"id2": "id1", "id1": "id2"}, inplace=True)

        merged_df = pd.concat([df, df_reverse], axis=0)
        merged_df.drop_duplicates(subset=["id1", "id2"] ,inplace=True)
        merged_df.reset_index(inplace=True, drop=True)
        merged_df.rename(columns={col: "ks"}, inplace=True)
        return merged_df

    def run_coll_dotplot(self, query_length, ref_length):
        dict2 = {"color": "strand"}

        coll_df = pd.read_csv(self.input_file, header=0, index_col=None, sep="\t", comment="#", low_memory=False)
        coll_df["queryChr"] = coll_df["queryChr"].astype(str)
        coll_df["refChr"] = coll_df["refChr"].astype(str)

        query_df, ref_df, coll_df = self.prefix_rm(query_length, ref_length, coll_df)
        blank_df = base.get_blank_chr(query_df, ref_df)

        query_chr = query_df['chr'].to_list()
        ref_chr = ref_df['chr'].to_list()
        coll_df = coll_df[coll_df['queryChr'].isin(query_chr)]
        coll_df = coll_df[coll_df['refChr'].isin(ref_chr)]

        coll_df['queryChr'] = pd.Categorical(coll_df['queryChr'], categories=query_chr, ordered=True)
        coll_df['refChr'] = pd.Categorical(coll_df['refChr'], categories=ref_chr, ordered=True)
        blank_df['queryChr'] = pd.Categorical(blank_df['queryChr'], categories=query_chr, ordered=True)
        blank_df['refChr'] = pd.Categorical(blank_df['refChr'], categories=ref_chr, ordered=True)

        custom_colors = {"+": "red", "-": "blue"}
        if self.type == 'order':
            dict1 = {"x": "queryId", "y": "refId"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=2, alpha=1) +
                    geom_blank(blank_df, show_legend=False) +
                    scale_color_manual(values=custom_colors) +
                    scale_x_continuous(expand=(0, 0)) +
                    scale_y_continuous(expand=(0, 0)))
        else:
            coll_df['queryStart'] = coll_df['queryStart'].apply(lambda x: x / 1000000)
            coll_df['referenceStart'] = coll_df['referenceStart'].apply(lambda x: x / 1000000)
            blank_df['queryStart'] = blank_df['queryStart'].apply(lambda x: x / 1000000)
            blank_df['referenceStart'] = blank_df['referenceStart'].apply(lambda x: x / 1000000)
            dict1 = {"x": "queryStart", "y": "referenceStart"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=2, alpha=1) +
                    geom_blank(blank_df, show_legend=False) +
                    scale_color_manual(values=custom_colors) +
                    scale_x_continuous(labels=self.major_formatter, expand=(0, 0)) +
                    scale_y_continuous(labels=self.major_formatter, expand=(0, 0)))

        plot1 = plot + theme_grey(base_size=50) + self.my_theme + labs(x=f'{self.query_name}', y=f'{self.ref_name}')
        plot1 = plot1 + guides(**{'color': guide_legend(override_aes={'size': 12, 'alpha': 1})})
        
        plot1.save(self.output_file_name, width=int(self.plotnine_figure_width), height=int(self.plotnine_figure_height), units="mm", limitsize=False)
        logger.info(f"Saving {self.plotnine_figure_width} x {self.plotnine_figure_height} mm image.")
        logger.info(f"Filename: {self.output_file_name}.")

    def run_coll_ks_dotplot(self, query_length, ref_length):
        dict2 = {"color": "ks"}
        ks_df = self.read_ks(self.ks, self.col)

        coll_df = pd.read_csv(self.input_file, header=0, index_col=None, sep="\t", comment="#", low_memory=False)
        coll_df["queryChr"] = coll_df["queryChr"].astype(str)
        coll_df["refChr"] = coll_df["refChr"].astype(str)

        query_df, ref_df, coll_df = self.prefix_rm(query_length, ref_length, coll_df)
        blank_df = base.get_blank_chr(query_df, ref_df)

        query_chr = query_df['chr'].to_list()
        ref_chr = ref_df['chr'].to_list()
        coll_df = pd.merge(coll_df, ks_df, how="left", left_on=["refGene", "queryGene"], right_on=["id1", "id2"])
        coll_df = coll_df[(coll_df['ks'] <= float(self.ks_area[1])) & (coll_df['ks'] >= float(self.ks_area[0]))]
        coll_df = coll_df[coll_df['queryChr'].isin(query_chr)]
        coll_df = coll_df[coll_df['refChr'].isin(ref_chr)]
        coll_df['queryChr'] = pd.Categorical(coll_df['queryChr'], categories=query_chr, ordered=True)
        coll_df['refChr'] = pd.Categorical(coll_df['refChr'], categories=ref_chr, ordered=True)
        blank_df['queryChr'] = pd.Categorical(blank_df['queryChr'], categories=query_chr, ordered=True)
        blank_df['refChr'] = pd.Categorical(blank_df['refChr'], categories=ref_chr, ordered=True)

        if self.type == 'order':
            dict1 = {"x": "queryId", "y": "refId"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=2, alpha=1) +
                    geom_blank(blank_df, show_legend=False) +
                    scale_color_cmap(cmap_name="gist_rainbow") +
                    scale_x_continuous(expand=(0, 0)) +
                    scale_y_continuous(expand=(0, 0)))
        else:
            coll_df['queryStart'] = coll_df['queryStart'].apply(lambda x: x / 1000000)
            coll_df['referenceStart'] = coll_df['referenceStart'].apply(lambda x: x / 1000000)
            blank_df['queryStart'] = blank_df['queryStart'].apply(lambda x: x / 1000000)
            blank_df['referenceStart'] = blank_df['referenceStart'].apply(lambda x: x / 1000000)
            dict1 = {"x": "queryStart", "y": "referenceStart"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=2, alpha=1) +
                    geom_blank(blank_df, show_legend=False) +
                    scale_color_cmap(cmap_name="gist_rainbow") +
                    scale_x_continuous(labels=self.major_formatter, expand=(0, 0)) +
                    scale_y_continuous(labels=self.major_formatter, expand=(0, 0)))

        plot1 = plot + theme_grey(base_size=50) + self.my_theme + labs(x=f'{self.query_name}', y=f'{self.ref_name}')
        if not self.output_file_name.endswith(".pdf"):
            plot1 = plot1 + theme(legend_position='none')
            plot1.save(self.output_file_name, width=int(self.plotnine_figure_width), height=int(self.plotnine_figure_height), units="mm", limitsize=False)
            img = Image.open(self.output_file_name)
            fig, ax = plt.subplots(figsize=(float(self.plotnine_figure_width) / 25.4, float(self.plotnine_figure_height) / 25.4))
            ax.imshow(img)
            ax.axis('off')

            norm = plt.Normalize(vmin=coll_df['ks'].min(), vmax=coll_df['ks'].max())
            sm = plt.cm.ScalarMappable(cmap='gist_rainbow', norm=norm)
            sm.set_array(coll_df['ks'])
            cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0, fraction=0.1)
            cbar.ax.tick_params(labelsize=40, length=5)
            cbar.set_label(r'${K_s}$', labelpad=30, fontsize=50, ha='center', va='center')
            fig.savefig(self.output_file_name, bbox_inches='tight')
        else:
            plot1.save(self.output_file_name, width=int(self.plotnine_figure_width),
                       height=int(self.plotnine_figure_height), units="mm", limitsize=False)
        logger.info(f"Saving {(float(self.plotnine_figure_width) / 25.4):.2f} x {(float(self.plotnine_figure_height) / 25.4):.2f} inch image.")
        logger.info(f"Filename: {self.output_file_name}.")

    def run_blast_dotplot(self, query_length, ref_length):
        dict2 = {"color": "strand"}

        coll_df = pd.read_csv(self.input_file, header=None, index_col=None, sep="\t", comment="#", low_memory=False)
        coll_df.columns = ["refGene", "refChr", "refId", "referenceStart", "referenceEnd", "refStrand",
                           "queryGene",	"queryChr", "queryId", "queryStart", "queryEnd", "queryStrand", "identity"]
        coll_df["strand"] = np.where(coll_df["refStrand"] == coll_df["queryStrand"], '+', '-')
        coll_df["queryChr"] = coll_df["queryChr"].astype(str)
        coll_df["refChr"] = coll_df["refChr"].astype(str)

        query_df, ref_df, coll_df = self.prefix_rm(query_length, ref_length, coll_df)
        blank_df = base.get_blank_chr(query_df, ref_df)

        query_chr = query_df['chr'].to_list()
        ref_chr = ref_df['chr'].to_list()
        coll_df = coll_df[coll_df['queryChr'].isin(query_chr)]
        coll_df = coll_df[coll_df['refChr'].isin(ref_chr)]
        coll_df['queryChr'] = pd.Categorical(coll_df['queryChr'], categories=query_chr, ordered=True)
        coll_df['refChr'] = pd.Categorical(coll_df['refChr'], categories=ref_chr, ordered=True)
        blank_df['queryChr'] = pd.Categorical(blank_df['queryChr'], categories=query_chr, ordered=True)
        blank_df['refChr'] = pd.Categorical(blank_df['refChr'], categories=ref_chr, ordered=True)

        custom_colors = {"+": "red", "-": "blue"}
        if self.type == "order":
            dict1 = {"x": "queryId", "y": "refId"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=1, alpha=0.5) +
                    geom_blank(blank_df, show_legend=False) +
                    scale_color_manual(values=custom_colors) +
                    scale_x_continuous(expand=(0, 0)) +
                    scale_y_continuous(expand=(0, 0)))
        else:
            coll_df['queryStart'] = coll_df['queryStart'].apply(lambda x: x / 1000000)
            coll_df['referenceStart'] = coll_df['referenceStart'].apply(lambda x: x / 1000000)
            blank_df['queryStart'] = blank_df['queryStart'].apply(lambda x: x / 1000000)
            blank_df['referenceStart'] = blank_df['referenceStart'].apply(lambda x: x / 1000000)
            dict1 = {"x": "queryStart", "y": "referenceStart"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=0.8, alpha=0.2) +
                    geom_blank(blank_df, show_legend=False) +
                    scale_color_manual(values=custom_colors) +
                    scale_x_continuous(labels=self.major_formatter, expand=(0, 0)) +
                    scale_y_continuous(labels=self.major_formatter, expand=(0, 0)))

        plot1 = plot + theme_grey(base_size=50) + self.my_theme + labs(x=f'{self.query_name}', y=f'{self.ref_name}')
        plot1 = plot1 + guides(color = guide_legend(override_aes={'size': 12, 'alpha': 1}))
        plot1.save(self.output_file_name, width=int(self.plotnine_figure_width), height=int(self.plotnine_figure_height), units="mm", limitsize=False)
        logger.info(f"Saving {self.plotnine_figure_width} x {self.plotnine_figure_height} mm image.")
        logger.info(f"Filename: {self.output_file_name}.")

    def run_blast_identity_dotplot(self, query_length, ref_length):
        dict2 = {"color": "identity"}

        table_df = pd.read_csv(self.input_file, header=None, index_col=None, sep="\t", comment="#", low_memory=False)
        table_df.columns = ["refGene", "refChr", "refId", "referenceStart", "referenceEnd", "refStrand",
                           "queryGene",	"queryChr", "queryId", "queryStart", "queryEnd", "queryStrand", "identity"]
        table_df["queryChr"] = table_df["queryChr"].astype(str)
        table_df["refChr"] = table_df["refChr"].astype(str)
        table_df["identity"] = table_df["identity"].astype(float) / 100

        query_df, ref_df, table_df = self.prefix_rm(query_length, ref_length, table_df)
        blank_df = base.get_blank_chr(query_df, ref_df)

        query_chr = query_df['chr'].to_list()
        ref_chr = ref_df['chr'].to_list()
        table_df = table_df[table_df['queryChr'].isin(query_chr)]
        table_df = table_df[table_df['refChr'].isin(ref_chr)]
        table_df['queryChr'] = pd.Categorical(table_df['queryChr'], categories=query_chr, ordered=True)
        table_df['refChr'] = pd.Categorical(table_df['refChr'], categories=ref_chr, ordered=True)
        blank_df['queryChr'] = pd.Categorical(blank_df['queryChr'], categories=query_chr, ordered=True)
        blank_df['refChr'] = pd.Categorical(blank_df['refChr'], categories=ref_chr, ordered=True)

        if self.type == "order":
            dict1 = {"x": "queryId", "y": "refId"}
            plot = (ggplot(table_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=1, alpha=0.5) +
                    scale_color_cmap(cmap_name="gist_rainbow") +
                    geom_blank(blank_df, show_legend=False) +
                    scale_x_continuous(expand=(0, 0)) +
                    scale_y_continuous(expand=(0, 0)))
        else:
            table_df['queryStart'] = table_df['queryStart'].apply(lambda x: x / 1000000)
            table_df['referenceStart'] = table_df['referenceStart'].apply(lambda x: x / 1000000)
            blank_df['queryStart'] = blank_df['queryStart'].apply(lambda x: x / 1000000)
            blank_df['referenceStart'] = blank_df['referenceStart'].apply(lambda x: x / 1000000)
            dict1 = {"x": "queryStart", "y": "referenceStart"}
            plot = (ggplot(table_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=0.8, alpha=0.2) +
                    scale_color_cmap(cmap_name="gist_rainbow") +
                    geom_blank(blank_df, show_legend=False) +
                    scale_x_continuous(labels=self.major_formatter, expand=(0, 0)) +
                    scale_y_continuous(labels=self.major_formatter, expand=(0, 0)))

        plot1 = plot + theme_grey(base_size=50) + self.my_theme + labs(x=f'{self.query_name}', y=f'{self.ref_name}')
        if not self.output_file_name.endswith(".pdf"):
            plot1 = plot1 + theme(legend_position='none')
            plot1.save(self.output_file_name, width=int(self.plotnine_figure_width), height=int(self.plotnine_figure_height), units="mm", limitsize=False)

            img = Image.open(self.output_file_name)
            fig, ax = plt.subplots(figsize=(float(self.plotnine_figure_width) / 25.4, float(self.plotnine_figure_height) / 25.4))
            ax.imshow(img)
            ax.axis('off')

            norm = plt.Normalize(vmin=table_df['identity'].min(), vmax=table_df['identity'].max())
            sm = plt.cm.ScalarMappable(cmap='gist_rainbow', norm=norm)
            sm.set_array(table_df['identity'])
            cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0, fraction=0.1)
            cbar.ax.tick_params(labelsize=40, length=5)
            cbar.set_label('identity', labelpad=30, fontsize=50, ha='center', va='center')
            fig.savefig(self.output_file_name, bbox_inches='tight')
            logger.info(f"Saving {(float(self.plotnine_figure_width) / 25.4):.2f} x {(float(self.plotnine_figure_height) / 25.4):.2f} inch image.")
            logger.info(f"Filename: {self.output_file_name}.")
        else:
            plot1.save(self.output_file_name, width=int(self.plotnine_figure_width), height=int(self.plotnine_figure_height), units="mm", limitsize=False)
        logger.info(f"Saving {(float(self.plotnine_figure_width) / 25.4):.2f} x {(float(self.plotnine_figure_height) / 25.4):.2f} inch image.")
        logger.info(f"Filename: {self.output_file_name}.")

    def run_coll_identity_dotplot(self, query_length, ref_length):
        dict2 = {"color": "score"}

        coll_df = pd.read_csv(self.input_file, header=0, index_col=None, sep="\t", comment="#", low_memory=False)
        coll_df["queryChr"] = coll_df["queryChr"].astype(str)
        coll_df["refChr"] = coll_df["refChr"].astype(str)

        query_df, ref_df, coll_df = self.prefix_rm(query_length, ref_length, coll_df)
        blank_df = base.get_blank_chr(query_df, ref_df)

        query_chr = query_df['chr'].to_list()
        ref_chr = ref_df['chr'].to_list()
        coll_df = coll_df[coll_df['queryChr'].isin(query_chr)]
        coll_df = coll_df[coll_df['refChr'].isin(ref_chr)]
        coll_df['queryChr'] = pd.Categorical(coll_df['queryChr'], categories=query_chr, ordered=True)
        coll_df['refChr'] = pd.Categorical(coll_df['refChr'], categories=ref_chr, ordered=True)
        blank_df['queryChr'] = pd.Categorical(blank_df['queryChr'], categories=query_chr, ordered=True)
        blank_df['refChr'] = pd.Categorical(blank_df['refChr'], categories=ref_chr, ordered=True)

        if self.type == 'order':
            dict1 = {"x": "queryId", "y": "refId"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=2, alpha=1) +
                    geom_blank(blank_df, show_legend=False) +
                    scale_color_cmap(cmap_name="gist_rainbow") +
                    scale_x_continuous(expand=(0, 0)) +
                    scale_y_continuous(expand=(0, 0)))
        else:
            coll_df['queryStart'] = coll_df['queryStart'].apply(lambda x: x / 1000000)
            coll_df['referenceStart'] = coll_df['referenceStart'].apply(lambda x: x / 1000000)
            blank_df['queryStart'] = blank_df['queryStart'].apply(lambda x: x / 1000000)
            blank_df['referenceStart'] = blank_df['referenceStart'].apply(lambda x: x / 1000000)
            dict1 = {"x": "queryStart", "y": "referenceStart"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=2, alpha=1) +
                    geom_blank(blank_df, show_legend=False) +
                    scale_color_cmap(cmap_name="gist_rainbow") +
                    scale_x_continuous(labels=self.major_formatter, expand=(0, 0)) +
                    scale_y_continuous(labels=self.major_formatter, expand=(0, 0)))

        plot1 = plot + theme_grey(base_size=50) + labs(x=f'{self.query_name}', y=f'{self.ref_name}')
        if not self.output_file_name.endswith(".pdf"):
            plot1 = plot1 +self.my_theme + theme(legend_position='none')
            plot1.save(self.output_file_name, width=int(self.plotnine_figure_width), height=int(self.plotnine_figure_height), units="mm", limitsize=False)
            img = Image.open(self.output_file_name)
            fig, ax = plt.subplots(figsize=(float(self.plotnine_figure_width) / 25.4, float(self.plotnine_figure_height) / 25.4))
            ax.imshow(img)
            ax.axis('off')

            norm = plt.Normalize(vmin=coll_df['score'].min(), vmax=coll_df['score'].max())
            sm = plt.cm.ScalarMappable(cmap='gist_rainbow', norm=norm)
            sm.set_array(coll_df['score'])
            cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0, fraction=0.1)
            cbar.ax.tick_params(labelsize=40, length=5)
            cbar.set_label(r'identity', labelpad=30, fontsize=50, ha='center', va='center')
            fig.savefig(self.output_file_name, bbox_inches='tight')
            logger.info(f"Saving {(float(self.plotnine_figure_width) / 25.4):.2f} x {(float(self.plotnine_figure_height) / 25.4):.2f} inch image.")
            logger.info(f"Filename: {self.output_file_name}.")
        else:
            plot1 = plot1 + self.my_theme
            plot1.save(self.output_file_name, width=int(self.plotnine_figure_width),
                       height=int(self.plotnine_figure_height), units="mm", limitsize=False)
        logger.info(f"Saving {(float(self.plotnine_figure_width) / 25.4):.2f} x {(float(self.plotnine_figure_height) / 25.4):.2f} inch image.")
        logger.info(f"Filename: {self.output_file_name}.")

    def dotplot_init(self):
        print()
        for key, value in vars(self).items():
            if key != "conf" and key not in ["my_theme", "col", "use_identity", "ks", "ks_area"]:
                print(key, "=", value)
            if key == "ks_area" and self.ks:
                print(key, "=", value)
            if key == "use_identity":
                if not self.ks:
                    print(key, "=", value)
                else:
                    pass
            if key == "ks" and self.ks:
                print(key, "=", value)
        print()

        if not self.query_length:
            logger.error("Please specify your query species chromosome length file")
            sys.exit(1)
        if not self.ref_length:
            logger.error("Please specify your reference species chromosome length file")
            sys.exit(1)
        if not self.input_file:
            logger.error("Please specify your input file(table file or collinearity file)")
            sys.exit(1)
        if not self.output_file_name:
            logger.error("Please specify your output file name")
            sys.exit(1)
        if self.disable_axis_text:
            self.my_theme += theme(axis_text=element_blank(),
                                   axis_ticks=element_blank())
        else:
            self.my_theme += theme(axis_ticks_major=element_line(size=4),
                                   axis_text_x=element_text(rotation=-60, hjust=0, vjust=1),
                                   axis_text=element_text(size=40))

        base.file_empty(self.input_file)
        base.output_file_parentdir_exist(self.output_file_name, self.overwrite)

    def run(self):
        logger.info("Dotplot module init and the following parameters are config information.")
        self.dotplot_init()
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Bitstream Vera Serif', 'Computer Modern Roman',
                                       'New Century Schoolbook', 'Century Schoolbook L', 'Utopia',
                                         'ITC Bookman', 'Bookman', 'Nimbus Roman No9 L', 'Times', 'Palatino', 'Charter', 'serif']
        plt.rcParams['mathtext.default'] = 'regular'

        with open(self.input_file, 'r') as file:
            first_line = file.readline()
        if first_line.startswith("#"):
            if self.ks:
                self.ks_area = self.ks_area.split(",")
                self.run_coll_ks_dotplot(self.query_length, self.ref_length)
            else:
                if self.use_identity:
                    self.run_coll_identity_dotplot(self.query_length, self.ref_length)
                else:
                    self.run_coll_dotplot(self.query_length, self.ref_length)
        else:
            if self.use_identity:
                self.run_blast_identity_dotplot(self.query_length, self.ref_length)
            else:
                self.run_blast_dotplot(self.query_length, self.ref_length)
