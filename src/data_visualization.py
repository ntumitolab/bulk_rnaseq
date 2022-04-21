import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from enum import Enum
import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, Isomap, MDS, SpectralEmbedding, LocallyLinearEmbedding
from sklearn.preprocessing import StandardScaler
import pandas as pd
from itertools import chain
from pathlib import Path

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
        
    def _get_group(self, group_name):
        gb = self.label_df.groupby(group_name)["sample"].apply(lambda x: list(x))
        groups = {f"{group_name}_{i}": gs for i, gs in gb.iteritems()}
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
                groups=self._get_group(groups)
                
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
                **kwargs):
        plt.figure(figsize=(10, 10))
        data = self._data[data_name]
        if groups:
            if isinstance(groups, str):
                groups=self._get_group(groups)
                
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
        x_label, y_label = method + "_1" + (f"({reducer.explained_variance_ratio_[0]:.2f})" if method=="PCA" else ""), method + "_2" + (f"({reducer.explained_variance_ratio_[1]:.2f})" if method=="PCA" else "")
        
        emb_data = pd.DataFrame(data=emb[:, :2], columns=[x_label, y_label], index=data.columns)
        emb_data["group"] = emb_data.index.to_series().map({vi : k for k, v in groups.items() for vi in v})
        g = sns.scatterplot(data=emb_data, x=x_label, y=y_label, hue="group", hue_order=group_orders, palette="bright")
        if file_name is not None:
            plt.savefig(file_name, dpi=300)
        
        
    def add_data_by_path(self, path, name):
        self._data[name] = pd.read_csv(path, sep='\t', index_col=0).dropna(how="all").fillna(0)
        
    def add_data(self, data, name):
        self._data[name] = data
        