from collections import namedtuple
from copy import copy, deepcopy
from typing import Any, Dict, List

import numpy as np
import tensorflow as tf
import umap

import seaborn as sns
import plotly.express as px
from visdom import Visdom

from xclone.model.XClone import XCloneVB


class FigureSpec:
    def __init__(self, figure_name: str, figure_kwargs: Dict[Any, Any] = None):
        if figure_kwargs is None:
            figure_kwargs = dict()
        self.fig_name = figure_name
        self.fig_kwargs = deepcopy(figure_kwargs)


class VisdomConvergenceTracker:
    """
    Domain-specific plots for real-time convergence diagnostics.
    Supported plot types: TODO
    """

    ####################################################################################################################
    #                                                 PUBLIC METHODS                                                   #
    ####################################################################################################################

    def __init__(self, visdom_obj: Visdom, xclone_instance: XCloneVB, figure_specs: List[FigureSpec] = None):
        self._visdom_obj = visdom_obj

        self._figure_specs = figure_specs
        self._validate_figure_specs()

        self._model = xclone_instance
        self._cnv_states = [
            f"({cnv[0]}, {cnv[1]})"
            for cnv in self._model.cnv_states.numpy().astype(np.int32)
        ]
        self._embedding = None
        self._precompute_embedding()

        self._cmap = "Viridis"
        self._cell_order = None
        self._deduce_cell_ordering()

        self._static_figures = {
            "count_amplification_heatmap",
            "reference_bafs_heatmap",
            "reference_total_cn_heatmap"
        }

    def ari_line_plot(self, ari_history):
        self._visdom_obj.line(
            Y=ari_history,
            X=np.arange(len(ari_history)),
            opts=dict(
                title="ARI",
                xlabel="iterations",
                ylabel="value",
                ymin=0,
                ymax=1,
                webgl=True
            ),
            win="ari"
        )

    def bafs_expected_from_cnv_heatmap(self):
        inferred_cnv = self._model.inferred_cnv
        inferred_clonal_labels = self._model.inferred_clonal_labels
        expected_baf_mx = []
        for clone in range(self._model.n_clones):
            clonal_cnv = inferred_cnv[:, clone]
            clonal_allelic_cn = tf.stack(tf.gather(self._model.cnv_states, clonal_cnv))
            clonal_total_cn = tf.reduce_sum(clonal_allelic_cn, axis=1)
            expected_clonal_bafs = clonal_allelic_cn[:, 0] / clonal_total_cn
            expected_baf_mx.append(expected_clonal_bafs)
        expected_baf_mx = np.column_stack(expected_baf_mx)
        sc_expected_baf_mx = np.column_stack([
            expected_baf_mx[:, cell_clone]
            for cell_clone in inferred_clonal_labels
        ])
        if self._cell_order is not None:
            sc_expected_baf_mx = sc_expected_baf_mx[:, self._cell_order]
        self._visdom_obj.heatmap(
            X=sc_expected_baf_mx,
            opts=dict(
                colormap=self._cmap,
                title="Bi-allele frequencies expected from CNV states",
                xlabel="cells",
                ylabel="segments",
            ),
            win="expected_ase"
        )

    def cell_projection_scatterplot(self):
        inferred_clone_labels = np.argmax(self._model.clonal_label_probs.mean().numpy(), axis=1)
        self._visdom_obj.scatter(
            X=self._embedding,
            Y=inferred_clone_labels + 1,
            opts=dict(
                title="UMAP projection colored by inferred clones",
                use_webgl=True
            ),
            win="proj"
        )

    def clone_probability_heatmap(self):
        inferred_clone_probs = self._model.clonal_label_probs.mean().numpy()
        self._visdom_obj.heatmap(
            X=inferred_clone_probs,
            opts=dict(
                colormap='Electric',
                title="Inferred clone probabilities",
                xlabel="clones",
                ylabel="cells",
                xmin=0,
                xmax=1
            ),
            win="probs"
        )

    def count_amplification_heatmap(self):
        counts_mx = self._model.total_counts_all_snps
        read_mapping_probs_mx = counts_mx / counts_mx.sum(axis=0)
        baseline = self._model.baseline_expression.numpy()
        amplification_mx = (1. / baseline).reshape(-1, 1) * read_mapping_probs_mx
        log2_fold_change_mx = np.log2(amplification_mx)
        log2_upper_bound = np.nanquantile(log2_fold_change_mx, .99)
        log2_lower_bound = np.nanquantile(log2_fold_change_mx, .01)
        overflow_mask = log2_fold_change_mx > log2_upper_bound
        log2_fold_change_mx[overflow_mask] = log2_upper_bound
        underflow_mask = log2_fold_change_mx < log2_lower_bound
        log2_fold_change_mx[underflow_mask] = log2_lower_bound

        self._visdom_obj.heatmap(
            X=log2_fold_change_mx[:, self._cell_order],
            opts=dict(
                colormap=self._cmap,
                title="log2(observed counts / baseline)",
                vmin=-3,
                vmax=3,
                xlabel="cells",
                ylabel="segments"
            ),
            win="count_amplification"
        )

    def count_amplification_vs_inferred_cnv_heatmap(self):
        raise NotImplementedError()

    def inferred_bafs_heatmap(self):
        inferred_cnv = self._model.inferred_cnv
        inferred_ase = self._model.xc_module_dict["ase"].bafs.mean().numpy()
        # print(inferred_ase)
        inferred_clonal_labels = np.argmax(self._model.clonal_label_probs.mean().numpy(), axis=1)
        inferred_baf_mx = []
        for clone in range(self._model.n_clones):
            clonal_cnv = inferred_cnv[:, clone]
            indices = clonal_cnv
            full_indices = tf.stack([tf.range(self._model.n_segments, dtype=indices.dtype), indices], axis=1)
            inferred_clonal_bafs = tf.gather_nd(inferred_ase, full_indices)
            inferred_baf_mx.append(inferred_clonal_bafs)
        inferred_baf_mx = np.column_stack(inferred_baf_mx)
        sc_inferred_baf_mx = np.column_stack([
            inferred_baf_mx[:, cell_clone]
            for cell_clone in inferred_clonal_labels
        ])
        if self._cell_order is not None:
            sc_inferred_baf_mx = sc_inferred_baf_mx[:, self._cell_order]
        self._visdom_obj.heatmap(
            X=sc_inferred_baf_mx,
            opts=dict(
                colormap=self._cmap,
                title="Inferred bi-allele frequencies",
                vmin=0,
                vmax=1,
                xlabel="cells",
                ylabel="segments",
            ),
            win="inferred_ase"
        )

    def inferred_bafs_vs_cnv_boxplot(self):
        inferred_ase = self._model.xc_module_dict["ase"].bafs.mean().numpy()
        self._visdom_obj.boxplot(
            X=inferred_ase,
            opts=dict(
                title="Inferred allelic ratios ",
                legend=self._cnv_states,
                xmin=0,
                xmax=1
            ),
            win="ase_boxplot"
        )

    def inferred_total_cn_heatmap(self):
        inferred_cnv_states = self._model.inferred_cnv.numpy()
        inferred_total_cn = np.array([
            [sum(self._model.cnv_states[cnv_ij]) for cnv_ij in cnv_i]
            for cnv_i in inferred_cnv_states
        ])
        inferred_clone_labels = np.argmax(self._model.clonal_label_probs.mean().numpy(), axis=1)
        # single-cell total copy number matrix
        sc_total_cn_mx = np.column_stack([
            inferred_total_cn[:, clone_label]
            for clone_label in inferred_clone_labels
        ])
        if self._cell_order is not None:
            sc_total_cn_mx = sc_total_cn_mx[:, self._cell_order]
        self._visdom_obj.heatmap(
            X=sc_total_cn_mx,
            opts=dict(
                title="Inferred total copy numbers",
                xlabel="cells",
                ylabel="segments",
                colormap=self._cmap
            ),
            win="inferred_total_cn"
        )

    def inferred_allelic_cnv_heatmap(self):
        inferred_cnv_states = self._model.inferred_cnv.numpy()
        self._visdom_obj.heatmap(
            X=inferred_cnv_states,
            opts=dict(
                colormap='Set2',
                title="Inferred CNV states",
                xlabel="cells",
                ylabel="segments"
            ),
            win="cnv"
        )

    def inferred_amplification_rates_scatterplot(self):
        self._visdom_obj.scatter(
            X=np.column_stack((
                np.arange(len(self._cnv_states)),
                np.exp(self._model.xc_module_dict["rdr"].amplification_factors.mean())
            )),
            Y=np.arange(len(self._cnv_states)) + 1,
            opts=dict(
                title="The mean of exp(amplification_factors)",
                xlabel="CNV states",
                ylabel="",
                legend=self._cnv_states,
            ),
            win="amplification_factors"
        )

    # def F_heatmap(self):
    #     F = tf.squeeze(self.model.F(), 0).numpy()
    #     self.visdom_obj.heatmap(
    #         X=F,
    #         opts=dict(
    #             colormap='Viridis',
    #             title="clonal_rdr_",
    #             xlabel="clones",
    #             ylabel="segments"
    #         ),
    #         win="clonal_rdr_"
    #     )

    def expected_loh_heatmap(self):
        inferred_cnv = self._model.inferred_cnv.numpy()
        loh_states = [
            i for i, cnv_state
            in enumerate(self._model.cnv_states)
            if 0 in cnv_state
        ]
        inferred_clone_labels = np.argmax(self._model.clonal_label_probs.mean().numpy(), axis=1)
        # single-cell total copy number matrix
        sc_total_cn_mx = np.column_stack([
            inferred_cnv[:, clone_label]
            for clone_label in inferred_clone_labels
        ])
        if self._cell_order is not None:
            sc_total_cn_mx = sc_total_cn_mx[:, self._cell_order]
        loh_mask = np.isin(sc_total_cn_mx, loh_states).astype(np.int32)
        self._visdom_obj.heatmap(
            X=loh_mask,
            opts=dict(
                colormap='Electric',
                title="Assumed loss of heterozygosity",
                xlabel="cells",
                ylabel="segments"
            ),
            win="loh_cnv"
        )

    def loss_history_lineplot(self, loss_history):
        self._visdom_obj.line(
            Y=loss_history,
            X=np.arange(len(loss_history)),
            opts=dict(
                title="ELBO loss",
                xlabel="iterations",
                ylabel="value",
                webgl=True
            ),
            win="elbo"
        )

    def reference_bafs_heatmap(self):
        self._visdom_obj.heatmap(
            X=(self._model.alt_counts / self._model.total_counts)[:, self._cell_order],
            opts=dict(
                title="Allelic rates from data: AD / DP",
                xlabel="cells",
                ylabel="segments",
                colormap=self._cmap
            ),
            win="ref_bafs"
        )

    def reference_total_cn_heatmap(self):
        self._visdom_obj.heatmap(
            X=self._model.reference_cn[:, self._cell_order],
            opts=dict(
                title="Reference copy numbers",
                xlabel="cells",
                ylabel="segments",
                colormap=self._cmap
            ),
            win="ref_total_cn"
        )

    def update_figures(self):
        for spec in self._figure_specs:
            if spec.fig_name in self._static_figures and self._model.epoch > 0:
                continue
            getattr(self, spec.fig_name)(**spec.fig_kwargs)

    ####################################################################################################################
    #                                                 PRIVATE METHODS                                                  #
    ####################################################################################################################

    def _deduce_cell_ordering(self):
        if self._model.reference_cn is None:
            self._cell_order = np.arange(self._model.n_cells)
        else:
            clustermap_graph = sns.clustermap(self._model.reference_cn, metric="cityblock", row_cluster=False)
            self._cell_order = clustermap_graph.dendrogram_col.reordered_ind

    def _precompute_embedding(self):
        projection_mapping = umap.UMAP(metric="manhattan")
        eps = 0.01
        projection_mapping.fit((self._model.alt_counts / (self._model.total_counts + eps)).T)
        self._embedding = projection_mapping.embedding_

    def _validate_figure_specs(self):
        object_methods = [method_name for method_name in dir(self)
                          if callable(getattr(self, method_name))]
        requested_figures = [spec.fig_name for spec in self._figure_specs]
        assert np.all(np.isin(requested_figures, object_methods)),\
            f"{np.setdiff1d(requested_figures, object_methods)} are not valid plots types.\n"\
            "Currently supported figures: TODO"
