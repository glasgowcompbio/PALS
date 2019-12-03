import warnings

import pandas as pd
from loguru import logger

from .pathway_analysis import PALS


def run_experiment(experiment_name, data_source, case, control, n_samples, significant_column, n_iter):
    res = {
        'experiment_name': experiment_name,
        'data_source': data_source,
        'case': case,
        'control': control,
        'n_samples': n_samples,
        'n_iter': n_iter,
        'significant_column': significant_column,
        'experiment_results': {}
    }

    # vary the number of (mzML) samples, run the pathway analysis methods for n_iter
    experiment_results = res['experiment_results']
    for n_sample in n_samples:
        experiment_results[n_sample] = [] # to track the length of the resampled results
        while len(experiment_results[n_sample]) < n_iter: # keep generating data until n_iter
                try:
                    # resample columns and generate a data source from it
                    i = len(experiment_results[n_sample])
                    logger.info('n_sample=%d iter=%d PALS experiment=%s case=%s control=%s' % (
                    n_sample, i, experiment_name, case, control))
                    ds_resampled = data_source.resample(case, control, n_sample)

                    # run PALS on the resampled data
                    pals = PALS(ds_resampled)
                    pathway_df = pals.get_pathway_df()
                    ora_df = pals.get_ora_df()

                    # store the results
                    item = {
                        # 'data': ds_resampled, # TOO BIG!!
                        'PALS': pathway_df,
                        'ORA': ora_df
                    }
                    experiment_results[n_sample].append(item)
                except UserWarning as e:
                    logger.warning('Failed to generate good data due to %s, will try again' % (str(e)))

    return res


def evaluate_performance(results, experiment_name, threshold, N):
    """
    Definition of precision and recall at N:
    Precision@N = (# of recommended items @N that are relevant) / (# of recommended items @N)
    Recall@N = (# of recommended items @N that are relevant) / (total # of relevant items)

    So here we define:
    - full results = top-N significant pathways from running the method on the complete dataframe
    - partial results = top-N significant pathways from running the method on the sampled dataframe having S samples
    - TP = in partial results and in full results
    - FP = in partial results and not in full results
    - FN = not in partial results and in full results

    From here, we can compute:
    Precision = TP/(TP+FP)
    Recall = TP/(TP+FN)
    """
    logger.debug('Generating PALS full results')
    res = results[experiment_name]
    significant_column = res['significant_column']
    experiment_results = res['experiment_results']
    ds = res['data_source']
    pals = PALS(ds)
    performances = []

    # generate PALS full results
    method = 'PALS'
    pals_full_df = pals.get_pathway_df()
    ora_full_df = pals.get_ora_df()
    pals_full = _select_significant_entries(pals_full_df, significant_column, N, threshold)
    ora_full = _select_significant_entries(ora_full_df, significant_column, N, threshold)

    # evaluate the partial results w.r.t to the full results
    logger.debug('Evaluating partial results')
    n_samples = list(experiment_results.keys())
    for n_sample in n_samples:
        iterations = range(len(experiment_results[n_sample]))
        for i in iterations:
            logger.debug('n_sample %d iteration %d' % (n_sample, i))
            item = experiment_results[n_sample][i]
            # for PALS
            method = 'PALS'
            full = pals_full
            partial_df = item[method]
            partial = _select_significant_entries(partial_df, significant_column, N, threshold)
            performances.append((method, n_sample, i) + _compute_prec_rec_f1(full, partial))
            # for ORA
            method = 'ORA'
            full = ora_full
            partial_df = item[method]
            partial = _select_significant_entries(partial_df, significant_column, N, threshold)
            performances.append((method, n_sample, i) + _compute_prec_rec_f1(full, partial))

    logger.debug('Done!')
    performance_df = pd.DataFrame(performances,
                                  columns=['method', 'n_sample', 'i', 'TP', 'FP', 'FN', 'precision', 'recall', 'F1'])
    return performance_df


def _select_significant_entries(pathway_df, significant_column, N, threshold):
    df = pathway_df.sort_values(significant_column, ascending=True)
    df = df.rename(columns={significant_column: 'p_value'})
    df = df[['pw_name', 'p_value']]
    df = df[df['p_value'] < threshold]
    return df[0:N]


def _compute_prec_rec_f1(full, partial):
    pathways_full = set(full.index.values)
    pathways_partial = set(partial.index.values)
    TP_items = pathways_full.intersection(pathways_partial)
    FP_items = pathways_partial - pathways_full
    FN_items = pathways_full - pathways_partial
    TP = len(TP_items)
    FP = len(FP_items)
    FN = len(FN_items)
    prec = 0
    rec = 0
    f1 = 0
    try:
        prec = TP / float(TP + FP)
        rec = TP / float(TP + FN)
        f1 = (2 * prec * rec) / (prec + rec)
    except:
        logger.warning('Something wrong!! TP=%d FP=%d FN=%d prec=%f rec=%f f1=%f' % (TP, FP, FN, prec, rec, f1))
    logger.debug('TP_items = %s' % TP_items)
    logger.debug('FP_items = %s' % FP_items)
    logger.debug('FN_items = %s' % FN_items)
    return TP, FP, FN, prec, rec, f1
