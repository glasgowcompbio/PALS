import copy
import warnings

import numpy as np
import pandas as pd
from loguru import logger
from scipy.stats import hypergeom
from scipy.stats import ttest_ind
from statsmodels.sandbox.stats.multicomp import multipletests

from .common import SIGNIFICANT_THRESHOLD, is_comparison_used, Method


class ORA(Method):

    def __init__(self, data_source, case=None, control=None):
        """
        Creates a ORA analysis
        :param data_source: a DataSource object
        """
        self.data_source = copy.deepcopy(data_source)
        self.case = case
        self.control = control

    ####################################################################################################################
    # public methods
    ####################################################################################################################

    def get_pathway_df(self, correct_multiple_tests=True, standardize=True):
        """
        Main method to perform over-representation (ORA) analysis
        :param: whether to log the initial data
        :return: a dataframe containing pathway analysis results from ORA
        """
        logger.debug('Calculating ORA')
        measurement_df = self.data_source.change_zero_peak_ints()
        if standardize:
            scaled_data = np.log(np.array(measurement_df))

            # Put the scaled data back into df for further use
            sample_names = measurement_df.columns
            measurement_df[sample_names] = scaled_data

        # For all of the pathways get all of the peak IDs
        assert len(self.data_source.dataset_pathways) > 0, 'No pathways found in the dataset'
        t_test_list = []
        pathways = self.data_source.dataset_pathways
        for pw in pathways:
            pathway_row_ids = self.data_source.dataset_pathways_to_row_ids[pw]
            pw_name = self.data_source.pathway_dict[pw]['display_name']
            path_params = [pw, pw_name]
            column_names = ['mapids', 'pw_name']
            for comp in self.data_source.comparisons:
                if not is_comparison_used(comp, self.case, self.control):
                    continue
                comparison_samples = self.data_source.get_comparison_samples(comp)
                condition_1 = comparison_samples[0]
                condition_2 = comparison_samples[1]

                # Perform hypergeometric test to assess the significance of finding formulae from dataset in the pathway
                # of interest. Parameters:
                # M = population size = the number of unique formulae found in database AND found in all the pathways
                # n = the number of success in the population = the number of differentially expressed entities in M
                # N = sample size = the number of unique formulae found in database AND found in the pathway of interest
                # k = the number of drawn success = the number of differentially expressed entities in N
                # P(X=k) = [(n C k) (M-n C N-k)] / (M c N)

                M = self.data_source.pathway_unique_ids_count
                significant_formulae = self._get_significant_formulae(condition_1, condition_2, measurement_df)

                n = len(self.data_source._get_pathway_unique_ids().intersection(significant_formulae))

                pw_ds_f = self.data_source.pathway_to_unique_ids_dict[pw]
                N = len(pw_ds_f)

                k = len(pw_ds_f.intersection(significant_formulae))

                # https://github.com/scipy/scipy/issues/7837
                sf = hypergeom.sf(k - 1, M, n, N)
                # logger.debug('pw=%s comp=%s M=%d n=%d N=%d k=%d sf=%f' % (pw, comp['name'], M, n, N, k, sf))

                # the combined p-value column is just the same as the p-value since there's nothing to combine
                item = (sf,)
                path_params.extend(item)

                # column names are computed is in the loop, but actually we only need it to be computed once
                col_name = comp['name'] + ' p-value'
                item = (col_name,)
                column_names.extend(item)

            t_test_list.append(path_params)

        t_test = pd.DataFrame(t_test_list, columns=column_names).set_index(['mapids'])
        t_test.index.name = 'mapids'

        # correct for multiple testing
        if correct_multiple_tests:
            logger.debug('Correcting for multiple t-tests')
        else:
            logger.debug('Not correcting for multiple t-tests')

        all_dfs = []
        for comp in self.data_source.comparisons:
            if not is_comparison_used(comp, self.case, self.control):
                continue

            # we use the combined p-value column name to store the p-values corrected after multiple testing
            col_name = comp['name'] + ' p-value'
            if self.data_source.database_name is not None:
                comb_col_name = '%s %s %s' % (self.data_source.database_name, comp['name'], 'comb_p')
            else:
                comb_col_name = '%s %s' % (comp['name'], 'comb_p')

            # copy the existing p-values
            pvalues = t_test[col_name].copy()
            if correct_multiple_tests:
                # check if any NaN, if yes exclude them
                keep = pd.notnull(pvalues)
                df = pvalues[keep]

                # perform multiple t-test corrections using FDR-BH
                reject, pvals_corrected, _, _ = multipletests(df.values, method='fdr_bh')

                # set the results back, and rename the column
                df.values[:] = pvals_corrected
            else:
                df = pvalues

            df = pd.DataFrame(df)
            df = df.rename(columns={col_name: comb_col_name})
            all_dfs.append(df)

        # combine all the results across all comparisons
        all_dfs = pd.concat(all_dfs, axis=1)

        # merge original df with the multiple-tests df
        t_test = t_test.merge(all_dfs, left_index=True, right_index=True, how='left')
        t_test_filled = t_test.fillna(1.0)

        mapids = t_test_filled.index.values.tolist()
        cov_df = self.data_source._calculate_coverage_df(mapids)
        coverage_df = cov_df.reindex(t_test_filled.index)  # make sure dfs are in same order before merging

        # Merge the two dfs together
        pathway_df = pd.merge(t_test_filled, coverage_df, left_index=True, right_index=True, how='outer')

        # del pathway_df.index.name
        pathway_df.rename_axis(None, inplace=True)

        return pathway_df

    ####################################################################################################################
    # private methods
    ####################################################################################################################

    def _get_significant_formulae(self, condition_1, condition_2, measurement_df, pathway_row_ids=None):
        if pathway_row_ids is None:
            pathway_row_ids = measurement_df.index.values
        c1 = measurement_df.loc[pathway_row_ids, condition_1].values
        c2 = measurement_df.loc[pathway_row_ids, condition_2].values
        with warnings.catch_warnings():
            # quietly ignore all warnings that come from failed t-tests if the values in the two conditions are the same
            warnings.filterwarnings('ignore')
            statistics, p_value = ttest_ind(c1, c2, axis=1)
        assert len(p_value) == len(pathway_row_ids)

        # TODO: vectorise this properly
        formula_detected_list = []
        for i in range(len(p_value)):
            row_id = pathway_row_ids[i]
            p = p_value[i]
            if p < SIGNIFICANT_THRESHOLD:
                peak_formulae = list(self.data_source.dataset_row_id_to_unique_ids[row_id])
                formula_detected_list.extend(peak_formulae)
        return set(formula_detected_list)
