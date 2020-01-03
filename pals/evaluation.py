import pandas as pd
from loguru import logger

from .common import SIGNIFICANT_THRESHOLD
from .pathway_analysis import PALS


def run_experiment(experiment_name, data_source, case, control, proportions, significant_column, n_iter,
                   plage_weight, hg_weight):
    res = {
        'experiment_name': experiment_name,
        'data_source': data_source,
        'case': case,
        'control': control,
        'proportions': proportions,
        'n_iter': n_iter,
        'significant_column': significant_column,
        'experiment_results': {},
        'plage_weight': plage_weight,
        'hg_weight': hg_weight
    }

    # vary the number of (mzML) samples, run the pathway analysis methods for n_iter
    experiment_results = res['experiment_results']
    n_samples = proportions * data_source.get_measurements().shape[0]
    n_samples = n_samples.astype(int)
    for i in range(len(proportions)):
        prop = proportions[i]
        n_sample = n_samples[i]
        experiment_results[prop] = []  # to track the length of the resampled results
        while len(experiment_results[prop]) < n_iter:  # keep generating data until n_iter
            try:
                # resample columns and generate a data source from it
                i = len(experiment_results[prop])
                logger.info('prop=%f n_sample=%d iter=%d PALS experiment=%s case=%s control=%s' % (
                    prop, n_sample, i, experiment_name, case, control))
                ds_resampled = data_source.resample(n_sample, case=case, control=control, axis=0)

                # run PALS on the resampled data
                pals = PALS(ds_resampled, plage_weight=plage_weight, hg_weight=hg_weight)
                pathway_df = pals.get_pathway_df()
                ora_df = pals.get_ora_df()

                # store the results
                item = {
                    # 'data': ds_resampled, # TOO BIG!!
                    'PALS': pathway_df,
                    'ORA': ora_df
                }
                experiment_results[prop].append(item)
            except UserWarning as e:
                logger.warning('Failed to generate good data due to %s, will try again' % (str(e)))

    return res


def evaluate_performance(results, experiment_name, N=None):
    """
    Definition of precision and recall at N:
    Precision@N = (# of recommended items @N that are relevant) / (# of recommended items @N)
    Recall@N = (# of recommended items @N that are relevant) / (total # of relevant items)
    If N is None, then all items above the significant threshold are considered.

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
    plage_weight = res['plage_weight']
    hg_weight = res['hg_weight']
    pals = PALS(ds, plage_weight=plage_weight, hg_weight=hg_weight)
    performances = []

    # generate PALS full results
    method = 'PALS'
    pals_full_df = pals.get_pathway_df()
    ora_full_df = pals.get_ora_df()
    pals_full = _select_significant_entries(pals_full_df, significant_column, SIGNIFICANT_THRESHOLD, N)
    ora_full = _select_significant_entries(ora_full_df, significant_column, SIGNIFICANT_THRESHOLD, N)

    # evaluate the partial results w.r.t to the full results
    logger.debug('Evaluating partial results')
    proportions = list(experiment_results.keys())
    for prop in proportions:
        iterations = range(len(experiment_results[prop]))
        for i in iterations:
            logger.debug('prop %d iteration %d' % (prop, i))
            item = experiment_results[prop][i]
            # for PALS
            method = 'PALS'
            full = pals_full
            partial_df = item[method]
            partial = _select_significant_entries(partial_df, significant_column, SIGNIFICANT_THRESHOLD, N)
            performances.append((method, prop, i) + _compute_prec_rec_f1(full, partial))
            # for ORA
            method = 'ORA'
            full = ora_full
            partial_df = item[method]
            partial = _select_significant_entries(partial_df, significant_column, SIGNIFICANT_THRESHOLD, N)
            performances.append((method, prop, i) + _compute_prec_rec_f1(full, partial))

    logger.debug('Done!')
    performance_df = pd.DataFrame(performances,
                                  columns=['method', 'proportion', 'i', 'TP', 'FP', 'FN', 'precision', 'recall', 'F1'])
    return performance_df


def _select_significant_entries(pathway_df, significant_column, threshold, N):
    df = pathway_df.sort_values(significant_column, ascending=True)
    df = df.rename(columns={significant_column: 'p_value'})
    try:
        df = df[['pw_name', 'p_value', 'sf', 'unq_pw_F', 'tot_ds_F', 'F_coverage']]
    except KeyError:
        df = df[['pw_name', 'p_value', 'unq_pw_F', 'tot_ds_F', 'F_coverage']]
    df = df[df['p_value'] < threshold]
    if N is not None:
        return df[0:N]
    else:
        return df


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
        pass
        # logger.warning('Something wrong!! TP=%d FP=%d FN=%d prec=%f rec=%f f1=%f' % (TP, FP, FN, prec, rec, f1))
    # logger.debug('TP_items = %s' % TP_items)
    # logger.debug('FP_items = %s' % FP_items)
    # logger.debug('FN_items = %s' % FN_items)
    return TP, FP, FN, prec, rec, f1
