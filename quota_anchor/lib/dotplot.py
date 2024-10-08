from plotnine import *
import pandas as pd
import numpy as np
from . import base

class Dotplot:
    def __init__(self, config_pra, parameter):
        self.overwrite = False
        self.col = "ks_NG86"
        self.type = "order"
        self.ref_name = "Reference species"
        self.query_name = "Query species"
        self.plotnine_figure_width=1500 
        self.plotnine_figure_height=1200
        for i in config_pra.sections():
            if i == 'dotplot':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)
        print()
        for key, value in vars(self).items():
            if key != "conf":
                print(key, "=", value)

        print()
        # please modify parameter
        self.my_theme = theme(
                panel_background=element_rect(fill='white'),
                panel_border=element_rect(fill=None, color="black", linewidth=2, linetype="solid"),
                axis_title=element_text(size=90),
                axis_ticks_major_x=element_line(size=4),
                axis_ticks_major_y=element_line(size=4),
                axis_ticks=element_blank(),
                legend_text=element_text(size=90, ha='center'),
                legend_title=element_text(size=70, ha='center'),
                axis_text_x=element_text(rotation=-60, hjust=0, vjust=1),
                # axis_text=element_blank(),
                axis_text=element_text(size=40),
                strip_text=element_text(size=50),
                panel_spacing=0,
                strip_background=element_blank()
            )

    @staticmethod
    def major_formatter(breaks):
        return [f'{int(b)}M' for b in breaks]

    @staticmethod
    def read_length(conf):
        df = pd.read_csv(conf, sep="\t", header=0, index_col=None)
        df['chr'] = df['chr'].astype(str)
        chr_list = df['chr']
        return chr_list


    @staticmethod
    def read_ks(file, col):
        df = pd.read_csv(file, sep="\t", header=0, index_col=None, low_memory=False)
        df["id1"] = df["id1"].astype(str)
        df["id2"] = df["id2"].astype(str)
        df = df.loc[:, ["id1", "id2", col]]
        df_reverse = df.loc[:, ["id2", "id1", col]]
        df_reverse.rename(columns={"id2": "id1", "id1": "id2"}, inplace=True)
        merged_df = pd.concat([df, df_reverse], axis=0)
        merged_df.drop_duplicates(subset=["id1", "id2"] ,inplace=True)
        merged_df.reset_index(inplace=True, drop=True)
        merged_df.rename(columns={col: "ks"}, inplace=True)
        return merged_df


    def run_coll_dotplot(self):
        dict2 = {"color": "strand"}
        query_chr = self.read_length(self.query_length)
        ref_chr = self.read_length(self.ref_length)
        coll_df = pd.read_csv(self.input_file, header=0, sep="\t", comment="#", low_memory=False)
        coll_df["queryChr"] = coll_df["queryChr"].astype(str)
        coll_df["refChr"] = coll_df["refChr"].astype(str)
        coll_df = coll_df[coll_df['queryChr'].isin(query_chr)]
        coll_df = coll_df[coll_df['refChr'].isin(ref_chr)]
        coll_df['queryChr'] = pd.Categorical(coll_df['queryChr'], categories=query_chr, ordered=True)
        coll_df['refChr'] = pd.Categorical(coll_df['refChr'], categories=ref_chr, ordered=True)
        blank_df = base.get_blank_chr(self.query_length, self.ref_length)
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
            
            plot1 = plot + theme_grey(base_size=50) + labs(x=f'{self.query_name}', y=f'{self.ref_name}')
            # please modify parameter
            plot1 = plot1 + guides(**{'color': guide_legend(override_aes={'size': 12, 'alpha': 1})}) + self.my_theme
            plot1.save(str(self.output_file_name), width=int(self.plotnine_figure_width), height=int(self.plotnine_figure_height), units="mm", limitsize=False)
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
            plot1.save(str(self.output_file_name), width=int(self.plotnine_figure_width), height=int(self.plotnine_figure_height), units="mm", limitsize=False)

    def run_coll_ks_dotplot(self):
        dict2 = {"color": "ks"}
        ks_df = self.read_ks(self.ks, self.col)
        query_chr = self.read_length(self.query_length)
        ref_chr = self.read_length(self.ref_length)
        coll_df = pd.read_csv(self.input_file, header=0, sep="\t", comment="#", low_memory=False)
        coll_df = pd.merge(coll_df, ks_df, how="left", left_on=["refGene", "queryGene"], right_on=["id1", "id2"])
        coll_df = coll_df[(coll_df['ks'] <= 4) & (coll_df['ks'] >= 0)]
        coll_df["queryChr"] = coll_df["queryChr"].astype(str)
        coll_df["refChr"] = coll_df["refChr"].astype(str)
        coll_df = coll_df[coll_df['queryChr'].isin(query_chr)]
        coll_df = coll_df[coll_df['refChr'].isin(ref_chr)]
        coll_df['queryChr'] = pd.Categorical(coll_df['queryChr'], categories=query_chr, ordered=True)
        coll_df['refChr'] = pd.Categorical(coll_df['refChr'], categories=ref_chr, ordered=True)
        blank_df = base.get_blank_chr(self.query_length, self.ref_length)

        if self.type == 'order':
            dict1 = {"x": "queryId", "y": "refId"}
            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=2, alpha=1) +
                    geom_blank(blank_df, show_legend=False) + 
                    scale_color_cmap(cmap_name="gist_rainbow") +
                    scale_x_continuous(expand=(0, 0)) +
                    scale_y_continuous(expand=(0, 0)))
            
            plot1 = plot + theme_grey(base_size=50) + labs(x=f'{self.query_name}', y=f'{self.ref_name}')
            # please modify parameter 
            plot1 = plot1 + guides(**{'color': guide_colorbar(display="gradient")}) +self.my_theme 
            plot1.save(str(self.output_file_name), width=int(self.plotnine_figure_width), height=int(self.plotnine_figure_height), units="mm", limitsize=False)
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
            plot1 = plot1 + guides(**{'color': guide_colorbar(display="gradient", alpha=1)}) +self.my_theme
            plot1.save(str(self.output_file_name), width=int(self.plotnine_figure_width), height=int(self.plotnine_figure_height), units="mm", limitsize=False)

    def run_blast_dotplot(self):
        dict2 = {"color": "strand"}
        query_chr = self.read_length(self.query_length)
        ref_chr = self.read_length(self.ref_length)
        coll_df = pd.read_csv(self.input_file, header=None, sep="\t", comment="#", low_memory=False)
        coll_df.columns = ["refGene", "refChr", "refId", "referenceStart", "referenceEnd", "refStrand",
                           "queryGene",	"queryChr", "queryId", "queryStart", "queryEnd", "queryStrand", "identity"]
        coll_df["strand"] = np.where(coll_df["refStrand"] == coll_df["queryStrand"], '+', '-')
        coll_df["queryChr"] = coll_df["queryChr"].astype(str)
        coll_df["refChr"] = coll_df["refChr"].astype(str)
        coll_df = coll_df[coll_df['queryChr'].isin(query_chr)]
        coll_df = coll_df[coll_df['refChr'].isin(ref_chr)]
        coll_df['queryChr'] = pd.Categorical(coll_df['queryChr'], categories=query_chr, ordered=True)
        coll_df['refChr'] = pd.Categorical(coll_df['refChr'], categories=ref_chr, ordered=True)
        blank_df = base.get_blank_chr(self.query_length, self.ref_length)
        custom_colors = {"+": "red", "-": "blue"}
        if self.type == "order":
            dict1 = {"x": "queryId", "y": "refId"}
            # coll_df['refId'] = coll_df['refId'].apply(lambda x: x / 1000)
            # coll_df['queryId'] = coll_df['queryId'].apply(lambda x: x / 1000)

            plot = (ggplot(coll_df, aes(**dict1)) +
                    facet_grid('refChr~queryChr', scales="free", space="free") +
                    geom_point(aes(**dict2), size=1, alpha=0.5) +
                    geom_blank(blank_df, show_legend=False) + 
                    scale_color_manual(values=custom_colors) +
                    scale_x_continuous(expand=(0, 0)) +
                    scale_y_continuous(expand=(0, 0)))

            plot1 = plot + theme_grey(base_size=50) + self.my_theme + labs(x=f'{self.query_name}', y=f'{self.ref_name}')
            plot1 = plot1 + guides(**{'color': guide_legend(override_aes={'size': 12, 'alpha': 1})})
            plot1.save(str(self.output_file_name), width=int(self.plotnine_figure_width), height=int(self.plotnine_figure_height), units="mm", limitsize=False)
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
            plot1 = plot1 + guides(**{'color': guide_legend(override_aes={'size': 12, 'alpha': 1})})
            plot1.save(str(self.output_file_name), width=int(self.plotnine_figure_width), height=int(self.plotnine_figure_height), units="mm", limitsize=False)

    def run(self):
        base.file_empty(self.input_file)
        base.output_file_parentdir_exist(self.output_file_name, self.overwrite)
        with open(self.input_file, 'r') as file:
            first_line = file.readline()
        if first_line.startswith("#"):
            if hasattr(self, "ks"):
                self.run_coll_ks_dotplot()
            else:
                self.run_coll_dotplot()
        else:
            self.run_blast_dotplot()
