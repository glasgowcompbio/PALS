import copy
import timeit
import warnings
from random import shuffle

import numpy as np
import pandas as pd
from loguru import logger
from scipy.stats import combine_pvalues
from scipy.stats import genextreme
from scipy.stats import hypergeom
from scipy.stats import ttest_ind

from .common import NUM_RESAMPLES, PLAGE_WEIGHT, HG_WEIGHT, is_comparison_used, Method


class PLAGE(Method):

    def __init__(self, data_source, num_resamples=NUM_RESAMPLES, plage_weight=PLAGE_WEIGHT, hg_weight=HG_WEIGHT,
                 case=None, control=None):
        """
        Creates a PALS analysis
        :param data_source: a DataSource object
        :param min_replace: replace a group with all zero values with this min intensity
        :param num_resamples: the number of times to resample p-values
        :param plage_weight: the weight for PLAGE (intensity) component when combining p-values
        :param hg_weight: the weight for hypergeometric component when combining p-values
        """
        self.data_source = copy.deepcopy(data_source)
        self.num_resamples = num_resamples

        # Add one to the expected number of pathway formulas for sf calculations - 100% gives a zero sf value and
        # subsequently effects all of the subsequent calculations
        self.plage_weight = plage_weight
        self.hg_weight = hg_weight

        self.case = case
        self.control = control

    ####################################################################################################################
    # public methods
    ####################################################################################################################

    def get_pathway_df(self, resample=True, standardize=True, streamlit_pbar=None):
        """
        Main method to perform pathway analysis
        :param resample: whether to perform resampling
        :param standardize: whether to standardize data
        :return: a dataframe containing pathway analysis results from PALS
        """
        activity_df = self.get_plage_activity_df(standardize)
        if resample:
            with warnings.catch_warnings():
                # FIXME: not sure if this is the best thing to do
                warnings.filterwarnings('ignore')
                plage_df = self.set_up_resample_plage_p_df(activity_df, streamlit_pbar=streamlit_pbar)
        else:
            plage_df = self.set_up_plage_p_df(activity_df)
        pathway_df = self.calculate_hg_values(plage_df)
        return pathway_df

    def get_plage_activity_df(self, standardize=True):
        """
        Performs data normalisation and computes the PLAGE activity dataframe
        :return: PLAGE activity dataframe
        """
        if standardize:
            # https://stackoverflow.com/questions/15933741/how-do-i-catch-a-numpy-warning-like-its-an-exception-not-just-for-testing
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                try:
                    measurement_df = self.data_source.standardize_intensity_df()
                except UserWarning as e:
                    # raise the exception if we encounter:
                    # "UserWarning: Numerical issues were encountered when scaling the data and might not be solved.
                    # The standard deviation of the data is probably very close to 0."
                    raise (e)
        else:
            measurement_df = self.data_source.get_measurements()

        # Standardizing, testing data
        mean = np.round(measurement_df.values.mean(axis=1))
        variance = np.round(measurement_df.values.std(axis=1))
        logger.debug("Mean values of the rows in the DF is %s" % str(mean))
        logger.debug("Variance in the rows of the DF is %s" % str(variance))

        activity_df = self._calculate_pathway_activity_df(measurement_df)
        return activity_df

    def set_up_resample_plage_p_df(self, activity_df, streamlit_pbar=None):
        """
        Obtains a PLAGE dataframe with resampling
        :param activity_df: a PLAGE activity dataframe
        :return: a dataframe with resampled pvalues
        """
        logger.debug("Calculating plage p-values with resampling")
        all_pvalues = [activity_df.index, activity_df['pw name']]
        column_names = ['pw_name']
        for comp in self.data_source.comparisons:
            if not is_comparison_used(comp, self.case, self.control):
                continue

            logger.debug('Comparison %s' % comp['name'])
            column_names.append(comp['name'] + ' p-value')
            null_max_tvalues = []
            null_min_tvalues = []
            comparison_samples = self.data_source.get_comparison_samples(comp)
            start = timeit.default_timer()
            for iteration in range(self.num_resamples):
                if iteration % 100 == 0:
                    logger.debug('Resampling %d/%d' % (iteration, self.num_resamples))
                if streamlit_pbar is not None:
                    progress = int((iteration + 1) / self.num_resamples * 100)
                    streamlit_pbar.progress(progress)
                condition_1, condition_2 = self._permute_two_lists(comparison_samples[0], comparison_samples[1])
                permutation_tvalues = self._calculate_t_values(activity_df, condition_1, condition_2)
                null_max_tvalues.append(max(permutation_tvalues))
                null_min_tvalues.append(min(permutation_tvalues))
            stop = timeit.default_timer()
            # logger.debug('Total time %d' % (stop - start))
            tvalues = self._calculate_t_values(activity_df, comparison_samples[0], comparison_samples[1])
            pvalues = self._compare_resamples(tvalues, null_max_tvalues, null_min_tvalues)
            all_pvalues.append(pvalues)

        t_test_list = map(list, zip(*all_pvalues))
        t_test = pd.DataFrame(t_test_list).set_index([0])
        t_test.columns = column_names
        t_test.index.name = 'mapids'
        t_test_filled = t_test.fillna(1.0)

        mapids = t_test_filled.index.values.tolist()
        cov_df = self.data_source._calculate_coverage_df(mapids)
        coverage_df = cov_df.reindex(t_test_filled.index)  # make sure dfs are in same order before merging

        # Merge the two dfs together
        pathway_df = pd.merge(t_test_filled, coverage_df, left_index=True, right_index=True, how='outer')
        return pathway_df

    def set_up_plage_p_df(self, activity_df):
        """
        Obtains a PLAGE dataframe (without resampling)
        :param activity_df: a PLAGE activity dataframe
        :return: a df containing mapids, pathway names and plage p-values for each comparison along with the
        Formula coverage for a dataset
        """
        logger.info("Calculating plage p-values")
        t_test_list = []
        for pathway, row in activity_df.iterrows():
            name = row[0]
            path_params = [pathway, name]
            column_names = ['pw_name']
            for comp in self.data_source.comparisons:
                if not is_comparison_used(comp, self.case, self.control):
                    continue

                comparison_samples = self.data_source.get_comparison_samples(comp)
                condition_1 = comparison_samples[0]
                condition_2 = comparison_samples[1]
                c1 = activity_df.loc[pathway, condition_1].values
                c2 = activity_df.loc[pathway, condition_2].values
                path_params.append(list(ttest_ind(c1, c2))[1])
                column_names.append(comp['name'] + ' p-value')
            t_test_list.append(path_params)

        t_test = pd.DataFrame(t_test_list).set_index([0])
        t_test.columns = column_names
        t_test.index.name = 'mapids'
        t_test_filled = t_test.fillna(1.0)

        mapids = t_test_filled.index.values.tolist()
        cov_df = self.data_source._calculate_coverage_df(mapids)
        coverage_df = cov_df.reindex(t_test_filled.index)  # make sure dfs are in same order before merging

        # Merge the two dfs together
        pathway_df = pd.merge(t_test_filled, coverage_df, left_index=True, right_index=True, how='outer')
        return pathway_df

    def calculate_hg_values(self, pathway_df):
        """
        Adds hypegeometric p-values to the pathway df
        :param pathway_df: a dataframe containing PLAGE scores and coverage for pathways
        :return: pathway_df with the hg scores and combined p-values added
        """
        # logger.debug("Calculating the hyper-geometric p-values")
        # Calculate the hg scores and create a temporary df to merge with the main df
        p_value_list = []
        mapids = pathway_df.index.values.tolist()

        for mp in mapids:
            M = self.data_source.pathway_unique_ids_count
            tot_pw_f = pathway_df.loc[mp]['unq_pw_F']
            n = tot_pw_f

            N = self.data_source.pathway_dataset_unique_ids_count
            formula_detected = pathway_df.loc[mp]['tot_ds_F']
            k = formula_detected

            # Add one to the expected number of pathway formulas for sf calculations - 100% gives a zero sf value and
            # subsequently effects all of the subsequent calculations
            # n = tot_pw_f + PW_F_OFFSET
            # sf = hypergeom.sf(k, M, n, N)

            sf = hypergeom.sf(k - 1, M, n, N)
            # logger.debug('pw=%s M=%d n=%d N=%d k=%d sf=%f' % (mp, M, n, N, k, sf))

            exp_value = hypergeom.mean(
                self.data_source.pathway_unique_ids_count,
                tot_pw_f,
                self.data_source.pathway_dataset_unique_ids_count).round(2)
            p_value_list.append([mp, sf, exp_value])

        p_value = pd.DataFrame(p_value_list)
        p_value_df = p_value.set_index([0])
        # Exp_F is the number of formula that would be expected in this pathway by chance
        p_value_df.columns = ['sf', 'exp_F']

        # merge the p_value df with the original pathways_df
        pathway_df_merge = pd.merge(pathway_df, p_value_df[['sf', 'exp_F']], left_index=True, right_index=True,
                                    how='outer')
        pathway_df_merge.F_coverage = pathway_df_merge.F_coverage.round(2)
        pathway_df_merge['Ex_Cov'] = ((pathway_df_merge['exp_F']) / pathway_df_merge['unq_pw_F']) * 100
        pathway_df_merge.Ex_Cov = pathway_df_merge.Ex_Cov.round(2)

        # logger.debug("Calculating the combined p-values")

        # Make a combined_p df to merge with the main df
        column_names = []
        for comp in self.data_source.comparisons:
            if not is_comparison_used(comp, self.case, self.control):
                continue
            if self.data_source.database_name is not None:
                col_name = '%s %s %s' % (self.data_source.database_name, comp['name'], 'comb_p')
            else:
                col_name = '%s %s' % (comp['name'], 'comb_p')
            column_names.append(col_name)

        combine_p_list = []
        for mp in mapids:
            combine_p_pathway = [mp]
            for comp in self.data_source.comparisons:
                if not is_comparison_used(comp, self.case, self.control):
                    continue

                p_value_colname = comp['name'] + ' p-value'
                p_value = pathway_df_merge.loc[mp][p_value_colname]
                sf = pathway_df_merge.loc[mp]['sf']
                com_p = combine_pvalues([p_value, sf], 'stouffer', [self.plage_weight, self.hg_weight])
                combine_p_pathway.append(com_p[1])
            combine_p_list.append(combine_p_pathway)

        comb_p = pd.DataFrame(combine_p_list)
        comb_p_df = comb_p.set_index([0])
        # Exp_F is the number of formula that would be expected in this pathway by chance
        comb_p_df.columns = column_names
        pathway_df_final = pd.merge(pathway_df_merge, comb_p_df[column_names], left_index=True, right_index=True,
                                    how='outer')
        return pathway_df_final

    ####################################################################################################################
    # private methods
    ####################################################################################################################

    def _calculate_pathway_activity_df(self, measurement_df):
        """
        Calculates the pathway activity DF given a dataframe of standardized intensities
        :param measurement_df: a standardized dataframe (of peak intensites)
        :return: a DF with Pathway names (rows) and the SVD activity levels for the samples (columns)
        """
        pathways = self.data_source.dataset_pathways

        # For all of the pathways get all of the peak IDs
        pathway_activities = []
        pw_names = []
        for pw in pathways:
            row_ids = self.data_source.dataset_pathways_to_row_ids[pw]
            pathway_data = measurement_df.loc[row_ids]  # DF selected from peak IDs.
            w, d, c = np.linalg.svd(np.array(pathway_data))

            pw_name = self.data_source.pathway_dict[pw]['display_name']
            pw_act_list = []
            pw_act_list.append(pw)
            pw_names.append(pw_name)
            pw_act_list.extend(list(c[0]))
            pathway_activities.append(pw_act_list)

        activity_df = pd.DataFrame(pathway_activities).set_index([0])
        activity_df.columns = measurement_df.columns
        activity_df.index.name = "Pathway ids"
        activity_df.insert(0, 'pw name', pw_names)
        return activity_df

    def _permute_two_lists(self, list1, list2):
        l = list(list1) + list(list2)
        shuffle(l)
        slist1 = l[0:len(list1)]
        slist2 = l[len(list1):len(l)]
        return slist1, slist2

    @staticmethod
    def _tstats(a):
        m = np.mean(a, axis=1)
        n = a.shape[1]
        var = np.sum(np.square(a - m[:, np.newaxis]), axis=1) / (n - 1)

        return m, n, var

    def _calculate_t_values(self, activity_df, condition_1, condition_2):

        c1 = activity_df.loc[:, condition_1].values
        c2 = activity_df.loc[:, condition_2].values

        m1, n1, var1 = PLAGE._tstats(c1)
        m2, n2, var2 = PLAGE._tstats(c2)

        se_total = np.sqrt(var1 / n1 + var2 / n2)
        se_total[se_total == 0.0] = np.nan
        tvalues = (m1 - m2) / se_total

        return tvalues

    def _compare_resamples(self, tvalues, null_max_tvalues, null_min_tvalues):
        pvalues = []
        maxparams = genextreme.fit(null_max_tvalues)
        minparams = genextreme.fit([-x for x in null_min_tvalues])
        for tvalue in tvalues:
            pvalue = genextreme.sf(tvalue, *maxparams) if tvalue >= 0 else genextreme.sf(-tvalue, *minparams)
            pvalues.append(pvalue)
        return pvalues
