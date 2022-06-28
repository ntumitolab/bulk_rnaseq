import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from enum import Enum
import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, Isomap, MDS, SpectralEmbedding, LocallyLinearEmbedding
from sklearn.preprocessing import StandardScaler
import colorcet as cc
import pandas as pd
from itertools import chain
from pathlib import Path
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from matplotlib.collections import PatchCollection

def transform_func(x):
    return np.log2(x + 1)


class DataPlotter:
    reducer_dic = {"PCA": PCA, "tSNE": TSNE, "UMAP": umap.UMAP}
    
    def __init__(self, label_path="../data/all_info/labels_mar.csv"):
        if Path(label_path).suffix == ".csv":
            self.label_df = pd.read_csv(label_path, index_col=0)
        elif Path(label_path).suffix == ".xlsx":
            self.label_df = pd.read_excel(label_path)
        self._data = {}
        
    def _get_group(self, group_name, selected_samples="all"):
        
        if isinstance(selected_samples, str) and selected_samples == "all":
            gb = self.label_df.groupby(group_name)["sample"].apply(lambda x: list(x))
        else:
            gb = self.label_df[self.label_df["sample"].isin(selected_samples)].groupby(group_name)["sample"].apply(lambda x: list(x))
        groups = {f"{group_name}: {i}": gs for i, gs in gb.iteritems()}
        return groups
    
    @staticmethod
    def _handle_group_colors(groups, group_orders, colors, samples):
        sample_groups = {vi : k for k, v in groups.items() for vi in v}
        color_dict = dict(zip(group_orders, colors[:len(group_orders)]))
        return [color_dict[sample_groups[c]] for c in samples]
        
        
    def plot_cluster(self,
                     data_name, 
                     groups=None,
                     transform_func=None,
                     group_orders=None,
                     selected_groups=None,
                     selected_rows=None,
                     filter_std=1,
                     figsize=(10, 20),
                     top_ele_ratio=0.1,
                     palette="bright",
                     cbar_label=None,
                     row_cluster=False,
                     row_dendrogram=True,
                     file_name=None,
                     **kwargs
                    ):
        col_colors, row_colors = None, None
        data = self._data[data_name]

        if groups:
            if isinstance(groups, str):
                groups=self._get_group(groups, selected_samples=data.columns)
                
            if selected_groups:
                groups = {k: v for k, v in groups.items() if k in selected_groups}
                data = data.loc[:, set(data.columns) & set(chain(*tuple(groups.values())))]
            
            if group_orders is None:
                group_orders = list(groups.keys())            
            colors = sns.color_palette(palette)
            col_colors = self._handle_group_colors(groups, group_orders, colors, data.columns)
        
        data = data[data.T.std() > 1]
        if selected_rows is not None:
            data = data.loc[list(set(selected_rows) & set(data.index)), :]
            
        if transform_func is not None:
            data = transform_func(data)
        
        g = sns.clustermap(data=data,
                           figsize=figsize,
                           row_cluster=row_cluster,
                           fmt=".3f",
                           col_colors=col_colors,
                           cbar_pos=(0, 1-top_ele_ratio+0.01, 0.05, top_ele_ratio-0.01),
                           dendrogram_ratio=top_ele_ratio,
                           colors_ratio=(0.02, 0.4 / figsize[1]),
                           cbar_kws={"label": cbar_label},
                           **kwargs)
        if file_name is not None:
            g.savefig(file_name, dpi=300)
    
    @staticmethod
    def _prepare_dr_x(x, standardize=True, transform_func=None):
        x = transform_func(x.copy()) if transform_func is not None else x.copy()
        x = StandardScaler().fit_transform(x.values.T) if standardize else x.values.T
        return x
    
        
    def plot_dr(self,
                data_name, 
                groups=None,
                group_orders=None,
                method="PCA", 
                selected_groups=None,
                selected_rows=None,
                standardize=True, 
                transform_func=None, 
                file_name=None,
                palette="bright",
                **kwargs):
        plt.figure(figsize=(10, 10))
        data = self._data[data_name]
        if groups:
            if isinstance(groups, str):
                groups=self._get_group(groups, selected_samples=data.columns)
                
            if selected_groups:
                groups = {k: v for k, v in groups.items() if k in selected_groups}
                data = data.loc[:, set(data.columns) & set(chain(*tuple(groups.values())))]
            
            if group_orders is None:
                group_orders = list(groups.keys())
        if selected_rows is not None:
            data = data.loc[list(set(selected_rows) & set(data.index)), :]
        
        x = self._prepare_dr_x(data, standardize, transform_func)
        n_components = min(x.shape[0], x.shape[1]) if method == "PCA" else 2
        kwargs.update({"n_components": n_components})
        # if method == "tSNE":
        #     kwargs.update({"learning_rate": "auto"})
        reducer = self.reducer_dic[method](**kwargs)
        emb = reducer.fit_transform(x)
        x_label, y_label = method + "_1" + (f"({reducer.explained_variance_ratio_[0]:.2f})" if method=="PCA" else ""), \
                            method + "_2" + (f"({reducer.explained_variance_ratio_[1]:.2f})" if method=="PCA" else "")
        
        emb_data = pd.DataFrame(data=emb[:, :2], columns=[x_label, y_label], index=data.columns)
        emb_data["group"] = emb_data.index.to_series().map({vi : k for k, v in groups.items() for vi in v})
        g = sns.scatterplot(data=emb_data, x=x_label, y=y_label, hue="group", hue_order=group_orders, palette=palette if not isinstance(palette, list) else palette[:len(emb_data["group"].unique())])
        if file_name is not None:
            plt.savefig(file_name, dpi=300)
        
        
    def add_data_by_path(self, path, name):
        self._data[name] = pd.read_csv(path, sep='\t', index_col=0).dropna(how="all").fillna(0)
        
    def add_data(self, data, name):
        self._data[name] = data
        
        
class Heatmap:
    def __init__(self, data, figsize):
        self.data = data
        self.figsize = figsize
        self.axes = []
        self.row_label_kw_dics = {}
        self._has_outter_row = False
        self.tickslabel_fs = 7
    
    @staticmethod
    def remove_border(ax):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
    
    def _add_row_color_util(self, ax, is_outter, row_labels, label, row_orders=None, rotation=0, palette="bright"):
        self.remove_border(ax)
        ax.set_xticks([])
        ax.set_xlabel(label, rotation=rotation)
        for label in ax.get_yticklabels():
            label.set_fontsize(self.tickslabel_fs)
        
        if not is_outter:
            ax.axes.yaxis.set_visible(False)

        uniq_row_labels = {u: [i for i, v in enumerate(row_labels) if v == u] for u in (set([uni for uni in row_labels]) if not row_orders else row_orders)}
        
        cols = sns.palettes.color_palette(palette, n_colors=8) + sns.palettes.color_palette("Set3", n_colors=8) + sns.color_palette(cc.glasbey, n_colors=12)
        pal = dict(zip(uniq_row_labels.keys(), cols[:len(uniq_row_labels)]))
        
        for i, (ul, uv) in enumerate(uniq_row_labels.items()):
            pc = PatchCollection([Rectangle((0, v), 1, 1) for v in uv], facecolors=pal[ul])
            ax.add_collection(pc)
        return pal
            
    def add_row_color(self, row_labels, label, row_orders=None, is_outter=False, rotation=45, palette="bright"):
        if is_outter and self._has_outter_row:
            raise ValueError("Already has an outter row")
        self._has_outter_row = is_outter
        self.row_label_kw_dics[label] = {"row_labels": row_labels, "row_orders": row_orders, "is_outter": is_outter, "rotation": rotation, "palette": palette}

    @staticmethod
    def _draw_legends(ax, pal_dic):
        has_legend = False
        y = 0.9
        last_legend = None
        for i, (label, p_dic) in enumerate(pal_dic.items()):
            h = [Patch(facecolor=color, label=lb) for lb, color in p_dic.items()]
            last_legend = ax.legend(handles= h, title=label, ncol=min(len(h), 3), bbox_to_anchor=(0.2 + 0.35 * i, y), loc="upper right")
            ax.add_artist(last_legend)
        ax.add_artist(last_legend)
        
    def show(self, ncols=10):
        fig = plt.figure(figsize=self.figsize)
        
        gs = fig.add_gridspec(nrows=5, ncols=ncols, left=0.05, right=0.95,
                              hspace=0.5, wspace=0.05)
        
        legend_ax = fig.add_subplot(gs[4, :])
        self.remove_border(legend_ax)
        data_ax = fig.add_subplot(gs[:4, -(ncols - len(self.row_label_kw_dics)):])
        sns.heatmap(self.data, ax=data_ax)
        data_ax.axes.yaxis.set_visible(False)
        for label in data_ax.get_xticklabels():
            label.set_fontsize(self.tickslabel_fs)
        
        legend_ax.set_xticks([])
        legend_ax.set_yticks([])
        
        ocp_cols = 1
        pal_dic = {}
        for label, kw in self.row_label_kw_dics.items():
            row_ax = fig.add_subplot(gs[:4, ocp_cols if not kw["is_outter"] else 0], sharey=data_ax)
            pal_dic[label] = self._add_row_color_util(row_ax, label=label, **kw)
            ocp_cols += (1 if not kw["is_outter"] else 0)
        self._draw_legends(legend_ax, pal_dic)
        