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
        'PALS': None,
        'ORA': None
    }

    # vary the number of (mzML) samples, run the pathway analysis methods for n_iter
    results = {}
    for n_sample in n_samples:
        results[n_sample] = []
        while len(results[n_sample]) < n_iter:

            # https://stackoverflow.com/questions/15933741/how-do-i-catch-a-numpy-warning-like-its-an-exception-not-just-for-testing
            with warnings.catch_warnings():
                warnings.filterwarnings('error')

                try:
                    # resample columns and generate a data source from it
                    i = len(results[n_sample])
                    logger.info('n_sample=%d iter=%d PALS experiment=%s case=%s control=%s' % (
                    n_sample, i, experiment_name, case, control))
                    ds_resampled = data_source.resample(case, control, n_sample)

                    # run PALS on the resampled data
                    pals = PALS(ds_resampled)
                    pathway_df = pals.get_pathway_df()

                    # store the results
                    item = {
                        'data': ds_resampled,
                        'result': pathway_df
                    }
                    results[n_sample].append(item)

                except Warning:
                    # to handle
                    # "UserWarning: Numerical issues were encountered when scaling the data and might not be solved.
                    # The standard deviation of the data is probably very close to 0."
                    logger.warning('Failed to generate good data, will try again')

    res['PALS'] = results
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
    ds = res['data_source']
    significant_column = res['significant_column']

    # generate full results
    pals = PALS(ds)
    pathway_df = pals.get_pathway_df()
    full = _select_significant_entries(pathway_df, significant_column, N, threshold)

    # evaluate PALS results
    logger.debug('Evaluating partial results')
    method = 'PALS'
    pals = res[method]
    performances = []
    n_samples = pals.keys()
    for n_sample in n_samples:
        for i in range(len(pals[n_sample])):
            item = pals[n_sample][i]
            data = item['data']
            df = item['result']
            partial = _select_significant_entries(df, significant_column, N, threshold)
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
    TP = len(pathways_full.intersection(pathways_partial))
    FP = len(pathways_partial - pathways_full)
    FN = len(pathways_full - pathways_partial)
    prec = 0
    rec = 0
    f1 = 0
    try:
        prec = TP / float(TP + FP)
        rec = TP / float(TP + FN)
        f1 = (2 * prec * rec) / (prec + rec)
    except:
        logger.warning('Something wrong!! TP=%d FP=%d FN=%d prec=%f rec=%f f1=%f' % (TP, FP, FN, prec, rec, f1))
    return TP, FP, FN, prec, rec, f1
