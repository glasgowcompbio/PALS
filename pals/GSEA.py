import copy

import numpy as np
import pandas as pd
from gseapy.gsea import GSEA
from loguru import logger

from .common import is_comparison_used, GSEA_RANKING_SNR, NUM_RESAMPLES, Method


class MSEA(GSEA):
    def load_data(self, cls_vec):
        # use all the peaks in the dataframe
        # TODO: we should add an improved peak -> formula filtering here
        exprs = self.data.copy()
        df = exprs.select_dtypes(include=[np.number])
        return df

    def run(self):
        try:
            super().run()
        except PermissionError:
            # ignore errors on windows:
            # PermissionError: [WinError 32] The process cannot access the file because it is being used by another process
            pass


class GSEA(Method):

    def __init__(self, data_source, num_resamples=NUM_RESAMPLES, method=GSEA_RANKING_SNR, case=None, control=None):
        """
        Creates a GSEA analysis
        This GSEA implementation is based on gseaPy.
        :param data_source: a DataSource object
        :param random_sets: the number of random permutation to do for GSEA
        :param pbar: whether to show progress bar or not
        """
        logger.debug('GSEA initialised with num_resamples=%d and ranking_method=%s' % (num_resamples, method))
        self.data_source = copy.deepcopy(data_source)
        self.num_resamples = num_resamples
        self.method = method
        self.case = case
        self.control = control

    ####################################################################################################################
    # public methods
    ####################################################################################################################

    def get_pathway_df(self):
        """
        Main method to perform GSEA/MSEA analysis
        :return: a dataframe containing pathway analysis results from GSEA
        """
        logger.debug('Calculating GSEA')
        measurement_df = self.data_source.change_zero_peak_ints()

        annot_df = self.data_source.get_annotations()
        joined = pd.merge(left=measurement_df, right=annot_df, left_index=True, right_index=True)
        joined = joined.set_index('entity_id')
        unique_ids = [self.data_source._get_unique_id(x) for x in joined.index.values]
        joined.index = unique_ids
        joined = joined.drop_duplicates(keep='first').sort_index()

        # gene_sets is a dict. key is pw name, values are a list of entries in that pathway
        gene_sets = {}
        assert len(self.data_source.dataset_pathways) > 0, 'No pathways found in the dataset'
        pathways = list(self.data_source.dataset_pathways)
        for pw in pathways:
            pathway_row_ids = self.data_source.dataset_pathways_to_row_ids[pw]
            pw_unique_ids = []
            for row_id in pathway_row_ids:
                pw_unique_ids.extend(self.data_source.dataset_row_id_to_unique_ids[row_id])
            pw_unique_ids = list(set(pw_unique_ids))
            gene_sets[pw] = pw_unique_ids

        # run GSEA for all comparisons
        all_dfs = []
        for comp in self.data_source.comparisons:
            if not is_comparison_used(comp, self.case, self.control):
                continue
            case = comp['case']
            control = comp['control']
            logger.debug('Running comparison case=%s control=%s' % (case, control))
            pheno_cols = set(self.data_source.get_experimental_design()['groups'][case])
            df_cols = measurement_df.columns.values

            # for each comparison, we need to create C (phenotype labels)
            # Loop over df_cols and store an indicator into C.
            # Entries in C is 1 if that column belongs to the case group, otherwise it's a 0
            C = []
            for col in df_cols:
                if col in pheno_cols:
                    C.append(1)
                else:
                    C.append(0)
            C = np.array(C)

            # actually runs GSEA here
            data = joined
            cls = C.tolist()
            outdir = None
            min_size = 1
            max_size = 1000
            permutation_num = self.num_resamples
            weighted_score_type = 1
            permutation_type = 'phenotype'
            method = self.method
            ascending = True
            processes = 1
            figsize = (6.5, 6)
            format = 'pdf',
            graph_num = 20
            no_plot = True
            seed = None
            verbose = False

            msea = MSEA(data, gene_sets, cls, outdir, min_size, max_size, permutation_num,
                        weighted_score_type, permutation_type, method, ascending, processes,
                        figsize, format, graph_num, no_plot, seed, verbose)
            msea.run()

            # convert GSEA results to dataframe
            df = msea.res2d
            df = df.reset_index()
            selected = df[['Term', 'pval', 'fdr', 'es']]
            selected = selected.rename(columns={'Term': 'mapids'}).set_index('mapids')

            col_name = comp['name'] + ' p-value'
            es_colname = comp['name'] + ' ES_score'
            if self.data_source.database_name is not None:
                comb_col_name = '%s %s %s' % (self.data_source.database_name, comp['name'], 'comb_p')
            else:
                comb_col_name = '%s %s' % (comp['name'], 'comb_p')

            pathway_df = selected.rename(columns={
                'pval': col_name,
                'es': es_colname,
                'fdr': comb_col_name
            })
            all_dfs.append(pathway_df)

        # combine all the results across all comparisons
        combined_df = pd.concat(all_dfs, axis=1, sort=False)
        combined_df.index.name = 'mapids'

        # create a dataframe of pathway mapids and names
        pw_name_df = []
        for map_id in pathways:
            pw_name = self.data_source.pathway_dict[map_id]['display_name']
            pw_name_df.append((map_id, pw_name))
        pw_name_df = pd.DataFrame(pw_name_df, columns=['mapids', 'pw_name']).set_index(['mapids'])
        combined_df = pw_name_df.merge(combined_df, left_index=True, right_index=True)

        # add formula coverage information
        mapids = combined_df.index.values.tolist()
        cov_df = self.data_source._calculate_coverage_df(mapids)
        coverage_df = cov_df.reindex(combined_df.index)  # make sure dfs are in same order before merging

        # Merge the two dfs together
        pathway_df = pd.merge(combined_df, coverage_df, left_index=True, right_index=True, how='outer')

        # del pathway_df.index.name
        pathway_df.rename_axis(None, inplace=True)

        return pathway_df
