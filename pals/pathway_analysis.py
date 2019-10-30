import timeit
from random import shuffle

import numpy as np
import pandas as pd
from loguru import logger
from scipy.stats import combine_pvalues
from scipy.stats import genextreme
from scipy.stats import hypergeom
from scipy.stats import ttest_ind
from sklearn import preprocessing

NUM_RESAMPLES = 1000


class PALS(object):

    # The constructor just takes in the analysis and defines the project
    def __init__(self, data_source, min_intensity=5000, num_resamples=NUM_RESAMPLES):
        self.data_source = data_source
        self.min_intensity = min_intensity  # replace a group with all zero values with this min intensity
        self.num_resamples = num_resamples

        # Add one to the expected number of pathway formulas for sf calculations - 100% gives a zero sf value and
        # subsequently effects all of the subsequent calculations
        self.PW_F_OFFSET = 1
        self.PLAGE_WEIGHT = 5
        self.HG_WEIGHT = 1

    ####################################################################################################################
    # public methods
    ####################################################################################################################

    """Main method to perform pathway analysis
    :returns: a dataframe containing pathway analysis results"""

    def get_pathway_df(self, resample=False):
        activity_df = self.get_plage_activity_df()
        if resample:
            plage_df = self.set_up_resample_plage_p_df(activity_df)
        else:
            plage_df = self.set_up_plage_p_df(activity_df)
        pathway_df = self.calculate_hg_values(plage_df)
        return pathway_df

    """ Method to get setup and return the pathway activity dataframe
    :returns: Dataframe of pathways (rows), samples (columns) and the calulated activity from SVD
    """

    def get_plage_activity_df(self):
        int_df = self._standardize_intensity_df(self.data_source.int_df)
        plage_activity_df = self._calculate_pathway_activity_df(int_df)
        return plage_activity_df

    """Obtains a plage dataframe with resampling"""

    def set_up_resample_plage_p_df(self, activity_df):
        logger.info("Calculating plage p-values with resampling")
        all_pvalues = [activity_df.index, activity_df['pw name']]
        column_names = ['pw_name']
        for comp in self.data_source.comparisons:
            logger.debug('Comparison %s' % comp['name'])
            column_names.append(comp['name'] + ' p-value')
            null_max_tvalues = []
            null_min_tvalues = []
            comparison_samples = self.data_source.get_comparison_samples(comp)
            start = timeit.default_timer()
            for iteration in range(self.num_resamples):
                if iteration % 100 == 0:
                    logger.debug('Resampling %d/%d' % (iteration, self.num_resamples))
                condition_1, condition_2 = self._permute_two_lists(comparison_samples[0], comparison_samples[1])
                permutation_tvalues = self._calculate_t_values(activity_df, condition_1, condition_2)
                null_max_tvalues.append(max(permutation_tvalues))
                null_min_tvalues.append(min(permutation_tvalues))
            stop = timeit.default_timer()
            logger.debug('Total time %d' % (stop - start))
            tvalues = self._calculate_t_values(activity_df, comparison_samples[0], comparison_samples[1])
            pvalues = self._compare_resamples(tvalues, null_max_tvalues, null_min_tvalues)
            all_pvalues.append(pvalues)

        t_test_list = map(list, zip(*all_pvalues))
        t_test = pd.DataFrame(t_test_list).set_index([0])
        t_test.columns = column_names
        t_test.index.name = 'mapids'
        t_test_filled = t_test.fillna(1.0)

        mapids = t_test_filled.index.values.tolist()
        cov_df = self.calculate_coverage_df(mapids)
        coverage_df = cov_df.reindex(t_test_filled.index)  # make sure dfs are in same order before merging

        # Merge the two dfs together
        pathway_df = pd.merge(t_test_filled, coverage_df, left_index=True, right_index=True, how='outer')
        return pathway_df

    """ Set up a df containing mapids, pathway names and plage p-values for
    each comparison along with the Formula coverage for a dataset.
    Param: An activity matrix from the PLAGE analysis
    Returns: A dataframe containing the above parameters,
    """

    def set_up_plage_p_df(self, activity_df):
        logger.info("Calculating plage p-values")
        t_test_list = []
        for pathway, row in activity_df.iterrows():
            name = row[0]
            path_params = [pathway, name]
            column_names = ['pw_name']
            for comp in self.comparisons:
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
        cov_df = self.calculate_coverage_df(mapids)
        coverage_df = cov_df.reindex(t_test_filled.index)  # make sure dfs are in same order before merging

        # Merge the two dfs together
        pathway_df = pd.merge(t_test_filled, coverage_df, left_index=True, right_index=True, how='outer')
        return pathway_df

    """ A method to add the hypegeometric p-values to the pathway df
        Params: The df containing the plage scores and coverage for a pathway
        Returns: The pathway_df with the hg scores and combined p-values added.
    """

    def calculate_hg_values(self, pathway_df):
        logger.info("Calculating the hyper-geometric p-values")
        # Calculate the hg scores and create a temporary df to merge with the main df
        p_value_list = []
        mapids = pathway_df.index.values.tolist()
        kegg_pw_f, ds_pw_f = self.data_source.get_unique_pw_f()

        for mp in mapids:
            tot_pw_f = pathway_df.loc[mp]['unq_pw_F']
            formula_detected = pathway_df.loc[mp]['tot_ds_F']
            sf = self.get_sf(formula_detected, kegg_pw_f, tot_pw_f + self.PW_F_OFFSET, ds_pw_f)
            exp_value = hypergeom.mean(kegg_pw_f, tot_pw_f, ds_pw_f).round(2)
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

        logger.info("Calculating the combined p-values")

        # Make a combined_p df to merge with the main df
        combine_p_list = []
        column_names = [comp['name'] + ' comb_p' for comp in self.data_source.comparisons]
        for mp in mapids:
            combine_p_pathway = [mp]
            for comp in self.data_source.comparisons:
                p_value_colname = comp['name'] + ' p-value'
                p_value = pathway_df_merge.loc[mp][p_value_colname]
                sf = pathway_df_merge.loc[mp]['sf']
                com_p = combine_pvalues([p_value, sf], 'stouffer', [self.PLAGE_WEIGHT, self.HG_WEIGHT])
                combine_p_pathway.append(com_p[1])
            combine_p_list.append(combine_p_pathway)

        comb_p = pd.DataFrame(combine_p_list)
        comb_p_df = comb_p.set_index([0])
        # Exp_F is the number of formula that would be expected in this pathway by chance
        comb_p_df.columns = column_names
        pathway_df_final = pd.merge(pathway_df_merge, comb_p_df[column_names], left_index=True, right_index=True,
                                    how='outer')
        return pathway_df_final

    """Writes dataframe as csv"""

    def write_df(self, title, df, suffix=''):
        filename = title + suffix + '.csv'
        df.to_csv(filename)
        logger.info("Your file has been written to: %s ", filename)

    ####################################################################################################################
    # private methods
    ####################################################################################################################

    """Method to set the zero values in a DF and standardize across the samples
    :param: int_df: Dataframe of peak intensites (raw)
    :returns: DF with zero intensities replaced and the values standardized
    """

    def _standardize_intensity_df(self, int_df):
        # Change the 0.00 intensities in the matrix to useable values
        int_df = self._change_zero_peak_ints(int_df)

        logger.debug("Scaling the data across the sample: zero mean and unit variance")

        # standardize the data across the samples (zero mean and unit variance))
        scaled_data = np.log(np.array(int_df))
        mean_std = np.mean(np.std(scaled_data, axis=1))
        scaled_data = preprocessing.scale(scaled_data, axis=1) * mean_std
        # Put the scaled data back into df for further use
        sample_names = int_df.columns
        int_df[sample_names] = scaled_data

        # Standardizing, testing data
        mean = np.round(int_df.values.mean(axis=1))
        variance = np.round(int_df.values.std(axis=1))
        logger.debug("Mean values of the rows in the DF is %s" % str(mean))
        logger.debug("Variance in the rows of the DF is %s" % str(variance))
        return int_df

    """ A method to change a 'zero' entries in a dataframe.
    If all intensities in a (factor) group are zero, a min value is set.
    If there are > 1 and < number in group zero intensities, then the average of the non_zeros entries is calculated
    and used. Assuming the PiMP mzXML file names are unique
    :param peak_int_df: A dataframe of peak intensities with peak ids (rows) and samples (columns)
    :returns: No return, modifies peak_int_df.
    """

    def _change_zero_peak_ints(self, peak_int_df):
        # Get the min_intensity value set for the analysis
        logger.debug("Setting the zero intensity values in the dataframe")
        # Replace 0.0 with NaN for easier operations ahead
        peak_int_df[peak_int_df == 0.0] = None
        for group_name, samples in self.data_source.groups.items():
            # If all zero in group then replace with minimum
            peak_int_df.loc[peak_int_df.loc[:, samples].isnull().all(axis=1), samples] = self.min_intensity

            # Replace any other zeros with mean of group
            subset_df = peak_int_df.loc[:, samples]
            peak_int_df.loc[:, samples] = subset_df.mask(subset_df.isnull(), subset_df.mean(axis=1), axis=0)
        return peak_int_df

    """ Method to calculate the pathway activity DF given a dataframe of standardized intensities
    :param int_df: Takes in a standardized dataframe (of peak intensites)
    :returns: A DF with Pathway names (rows) and the SVD activity levels for the samples (columns)
    """

    def _calculate_pathway_activity_df(self, int_df):
        pathways = self.data_source.ds_pathways

        # For all of the pathways get all of the peak IDs
        pathway_activities = []
        pw_names = []
        for pw in pathways:
            peak_ids = self.data_source.ds_pathways_peak_ids[pw]
            pathway_peaks = int_df.loc[peak_ids] / np.sqrt(len(peak_ids))  # DF selected from peak IDs.
            w, d, c = np.linalg.svd(np.array(pathway_peaks))

            pw_name = self.data_source.pathway_dict[pw]
            pw_act_list = []
            pw_act_list.append(pw)
            pw_names.append(pw_name)
            pw_act_list.extend(list(c[0] * d[0]))

            pathway_activities.append(pw_act_list)

            activity_df = pd.DataFrame(pathway_activities).set_index([0])
            activity_df.columns = int_df.columns
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

        m1, n1, var1 = PALS._tstats(c1)
        m2, n2, var2 = PALS._tstats(c2)

        se_total = np.sqrt(var1 / n1 + var2 / n2)
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

    """ Calculate the Formula coverage for a dataset.
        Param: The Mapids for the pathways to be used
        Returns: A dataframe containing the number of unique formulae for a pathway, along with those
        annotated, identified and the total unique Fomulae in each pathway.
    """

    def calculate_coverage_df(self, mapids):
        logger.info("Calculating dataset formula coverage")

        # num_formula: Stores the number of unqique kegg formulae for a pathway
        num_formula = self.data_source.get_pw_unique_F(mapids)
        num_totalF = self.data_source.get_ds_pw_compounds(mapids)

        # unq_pw_f: unique formula expected for a pathway
        # tot_ds_F: unique formula for a pathway in a dataset
        data = {'unq_pw_F': num_formula,
                'tot_ds_F': num_totalF, }
        for_df = pd.DataFrame(data)
        for_df.index = mapids

        # Calculate the coverage of the formula found in the ds vs formulae in the pathway
        for_df['F_coverage'] = (((for_df['tot_ds_F']) / for_df['unq_pw_F']) * 100).round(2)

        return for_df

    # Get the Hypergeometric p-value
    def get_sf(self, formula_dectected, kegg_pw_f, tot_pw_f, ds_pw_f):
        [k, M, n, N] = [formula_dectected, kegg_pw_f, tot_pw_f, ds_pw_f]
        prb = hypergeom.sf(k, M, n, N)
        # print "This would suggest that there is a ", prb, "chance that ",formula_detected, "or greater metabolites would be detected randomly"
        return prb
