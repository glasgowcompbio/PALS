import numpy as np
import pandas as pd
from loguru import logger
from tqdm import tqdm

from .common import is_comparison_used


class GSEA(object):

    def __init__(self, data_source, random_sets=1000, p_exp=1, pbar=True, case=None, control=None):
        """
        Creates a GSEA analysis
        This GSEA implementation is based on https://github.com/mrcinv/GSEA.py.
        :param data_source: a DataSource object
        :param random_sets: the number of random permutation to do for GSEA
        :param p_exp: an exponent p to control the weight of the step in GSEA
        :param pbar: whether to show progress bar or not
        """
        self.data_source = data_source
        self.random_sets = random_sets
        self.pbar = pbar
        self.p_exp = p_exp

        self.case = case
        self.control = control


    ####################################################################################################################
    # public methods
    ####################################################################################################################

    def get_pathway_df(self, standardize=True):
        """
        Main method to perform GSEA/MSEA analysis
        :param: whether to log the initial data
        :return: a dataframe containing pathway analysis results from GSEA
        """
        logger.debug('Calculating GSEA')
        if standardize:
            measurement_df = self.data_source.change_zero_peak_ints()
            scaled_data = np.log(np.array(measurement_df))

            # Put the scaled data back into df for further use
            sample_names = measurement_df.columns
            measurement_df[sample_names] = scaled_data

        D = measurement_df.values

        # gene_sets is a list of lists.
        # Each list in gene_sets contains indices of rows in measurement_df for that pathway.
        gene_sets = []
        assert len(self.data_source.dataset_pathways) > 0, 'No pathways found in the dataset'
        pathways = list(self.data_source.dataset_pathways)
        for pw in pathways:
            pathway_row_ids = self.data_source.dataset_pathways_to_row_ids[pw]
            row_locs = list(
                set([measurement_df.index.get_loc(row_id) for row_id in pathway_row_ids]))  # remove duplicate?
            gene_sets.append(row_locs)
            # logger.debug('%s %s' % (pw, row_locs))

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
            order, NES, p_values = self.gsea(D, C, gene_sets, p_exp=self.p_exp, random_sets=self.random_sets,
                                             pbar=self.pbar)

            # convert GSEA results to dataframe
            data = []
            for idx, nes, pval in zip(order, NES, p_values):
                map_id = pathways[idx]
                pw_name = self.data_source.pathway_dict[map_id]['display_name']
                data.append([map_id, nes, pval, pval])

            # here the combined p-value is just the same as the p-value
            col_name = comp['name'] + ' p-value'
            es_colname = comp['name'] + ' ES_score'
            if self.data_source.database_name is not None:
                comb_col_name = '%s %s %s' % (self.data_source.database_name, comp['name'], 'comb_p')
            else:
                comb_col_name = '%s %s' % (comp['name'], 'comb_p')

            pathway_df = pd.DataFrame(data, columns=['mapids', es_colname, col_name, comb_col_name]).set_index(
                ['mapids'])
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
        del pathway_df.index.name
        return pathway_df

    def gsea(self, D, C, S_sets, p_exp=1, random_sets=1000, pbar=True):
        """Performs Multiple Hypotesis Testing.

        Arguments:
        ----------
        D: a 2D array of expression data for N genes and k sample

        C: a list of values 1, if a i-th sample is in the phenotype or 0 otherwise

        S_sets: a list of variable length of gene indexes that belong to a gene set

        p_exp: exponent parameter to control the weight of the step (defaults to 1)

        random_sets: number of randomly generated gene sets


        Returns:
        --------

        order: list of gene indexes, ordered by the greatest
            Normalized Enrichment Scores

        NES: a list of Normalized Enrichment Scores (ordered by absolute value)

        p_value: a list of p-values for the NES
        """
        N, k = D.shape
        n = len(S_sets)
        # generate random gene sets
        p_value = np.zeros(n)
        NES = np.zeros(n)
        ES = np.zeros(n)
        ES_pi = np.zeros((random_sets, n))
        L, r = self.rank_genes(D, C)
        # enrichment scores for S_i
        for i in range(n):
            ES[i] = self.enrichment_score(L, r, S_sets[i], p_exp)

        # whether to use progress bar
        iterable = range(random_sets)
        if pbar:
            iterable = tqdm(iterable)

        for i in iterable:
            pi = np.array([np.random.randint(0, 2) for i in range(k)])
            L, r = self.rank_genes(D, pi)
            ES_pi[i, :] = [self.enrichment_score(L, r, S_sets[j], p_exp)
                           for j in range(n)]

        # calculate normalized enrichment scores and p-values
        for i in range(n):
            # normalize separately positive and negative values
            ES_plus = ES_pi[:, i][ES_pi[:, i] > 0]
            ES_minus = ES_pi[:, i][ES_pi[:, i] < 0]
            mean_plus = np.mean(ES_plus)
            mean_minus = np.mean(ES_minus)
            if ES[i] > 0:
                NES[i] = ES[i] / mean_plus
                p_value[i] = sum(ES_plus > ES[i]) / len(ES_plus)
            elif ES[i] < 0:
                NES[i] = -ES[i] / mean_minus
                p_value[i] = sum(ES_minus < ES[i]) / len(ES_minus)

        NES_sort = sorted(enumerate(NES), key=lambda x: -abs(x[1]))
        order = [x[0] for x in NES_sort]
        NES = [x[1] for x in NES_sort]
        return order, NES, p_value[order]

    def rank_genes(self, D, C):
        """Ranks genes in expression dataset according to the correlation with a
        phenotype.

        Arguments:
        ----------
        D: a 2D array of expression data for N genes and k samples

        C: a list of values 1, if a i-th sample is in the phenotype or 0 otherwise

        Returns:
        --------
        L: an ordered list of gene indexes

        r: a ordered list of correlation coefficients

        """
        N, k = D.shape
        # way faster than np.corrcoef
        C = np.array(C)
        ED = np.mean(D, 1)
        EC = np.mean(C)
        EDC = np.mean(D * C, 1)
        KOV = EDC - ED * EC
        sD = (np.mean(D ** 2, 1) - ED ** 2) ** 0.5
        sC = (np.mean(C ** 2) - EC ** 2) ** 0.5
        rL = KOV / sD / sC
        # rL = []
        # for i in range(N):
        #     rL.append(np.corrcoef(D[i,:],C)[0,1])

        rL = sorted(enumerate(rL), key=lambda x: -x[1])
        r = [x[1] for x in rL]
        L = [x[0] for x in rL]
        return L, r

    def enrichment_score(self, L, r, S, p_exp):
        """Calculates enrichment score (ES) for a given gene expression data.

        Arguments:
        ---------
        L: an ordered list of indexes for N genes

        r: a list of correlations of a gene expression data with phenotypes for
            N genes

        C: a list of classes(phenotypes) for k samples

        S: a list of gene indexes that belong to a gene set

        p_exp: exponent parameter to control the weight of the step (defaults to 1)

        Returns:
        -------

        ES: enrichment score for the given gene set S
        """

        N = len(L)
        S_mask = np.zeros(N)
        S_mask[S] = 1
        # reorder gene set mask
        S_mask = S_mask[L]
        N_R = sum(abs(r * S_mask) ** p_exp)
        if N_R != 0:
            P_hit = np.cumsum(abs(r * S_mask) ** p_exp) / N_R
        else:
            P_hit = np.zeros_like(S_mask)
        N_H = len(S)
        if N != N_H:
            P_mis = np.cumsum((1 - S_mask)) / (N - N_H)
        else:
            P_mis = np.zeros_like(S_mask)
        idx = np.argmax(abs(P_hit - P_mis))
        return P_hit[idx] - P_mis[idx]
