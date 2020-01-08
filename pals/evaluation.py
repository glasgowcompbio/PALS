import warnings

import pandas as pd
from loguru import logger

from .GSEA import GSEA
from .ORA import ORA
from .PALS import PALS
from .common import SIGNIFICANT_THRESHOLD, MIN_REPLACE, set_log_level_info, set_log_level_debug, NUM_RESAMPLES, \
    PLAGE_WEIGHT, HG_WEIGHT
from .noise import construct_intensity_df, add_random_peaks, convert_to_data_source


########################################################################################################################
# Experiments with adding noise on synthetic data
########################################################################################################################

def run_noise_experiment(background_pathways, case_fnames, control_fnames, pathway_names, num_iterations, plage_weight,
                         hg_weight, gsea_resamples, reqd_scenarios, pbar=False, parallel=False):
    pals_dfs = []
    ora_dfs = []
    gsea_dfs = []
    random = False  # whether to sample intensity data randomly or using pre-set values

    for i in range(len(reqd_scenarios)):
        scenario = reqd_scenarios[i]
        logger.info(scenario)
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', r'divide by zero')
            warnings.filterwarnings('ignore', r'invalid value encountered')
            results = calc_av_p_scores(case_fnames, control_fnames, pathway_names, num_iterations,
                                       percent=scenario['percent'], random=random, noise_std=scenario['noise_std'],
                                       prob_missing_peaks=scenario['prob_missing_peaks'],
                                       background_pathways=background_pathways, plage_weight=plage_weight,
                                       hg_weight=hg_weight, gsea_resamples=gsea_resamples, pbar=pbar, parallel=parallel)

        df = construct_single_box_df(results, scenario['percent'], scenario['prob_missing_peaks'],
                                     scenario['noise_std'], 'PALS')
        pals_dfs.append(df)

        df = construct_single_box_df(results, scenario['percent'], scenario['prob_missing_peaks'],
                                     scenario['noise_std'], 'ORA')
        ora_dfs.append(df)

        df = construct_single_box_df(results, scenario['percent'], scenario['prob_missing_peaks'],
                                     scenario['noise_std'], 'GSEA')
        gsea_dfs.append(df)

    pals_df = pd.concat(pals_dfs, axis=0)
    ora_df = pd.concat(ora_dfs, axis=0)
    gsea_df = pd.concat(gsea_dfs, axis=0)

    return pals_df, ora_df, gsea_df


def calc_av_p_scores(case_fnames, control_fnames, pathway_names, num_iterations, percent=0, random=False, noise_mean=0,
                     noise_std=5, prob_missing_peaks=0.2, background_pathways=100, min_replace=MIN_REPLACE,
                     plage_weight=PLAGE_WEIGHT, hg_weight=HG_WEIGHT, gsea_resamples=NUM_RESAMPLES, parallel=False,
                     pbar=False):
    sample_fnames = control_fnames + case_fnames
    set_log_level_info()

    # create num_iterations random data
    params = _get_noise_params(background_pathways, case_fnames, control_fnames, gsea_resamples, hg_weight,
                               min_replace, noise_mean, noise_std, pathway_names, pbar, percent, plage_weight,
                               prob_missing_peaks, random, sample_fnames)
    if not parallel:
        noise_results = []
        for it in range(num_iterations):
            item = _get_noise_performance(params)
            noise_results.append(item)
    else:
        all_params = [params for i in range(num_iterations)]
        import ipyparallel as ipp
        rc = ipp.Client()
        dview = rc[:]  # use all engines​
        with dview.sync_imports():
            pass
        noise_results = dview.map_sync(_get_noise_performance, all_params)

    # returns the final results
    results = {
        'PALS': [item['PALS'] for item in noise_results],
        'ORA': [item['ORA'] for item in noise_results],
        'GSEA': [item['GSEA'] for item in noise_results]
    }
    set_log_level_debug()
    return results


def _get_noise_params(background_pathways, case_fnames, control_fnames, gsea_resamples, hg_weight,
                      min_replace, noise_mean, noise_std, pathway_names, pbar, percent, plage_weight,
                      prob_missing_peaks, random, sample_fnames):
    return {
        'background_pathways': background_pathways,
        'case_fnames': case_fnames,
        'control_fnames': control_fnames,
        'gsea_resamples': gsea_resamples,
        'hg_weight': hg_weight,
        'min_replace': min_replace,
        'noise_mean': noise_mean,
        'noise_std': noise_std,
        'pathway_names': pathway_names,
        'pbar': pbar,
        'percent': percent,
        'plage_weight': plage_weight,
        'prob_missing_peaks': prob_missing_peaks,
        'random': random,
        'sample_fnames': sample_fnames
    }


def _get_noise_performance(params):
    background_pathways = params['background_pathways']
    case_fnames = params['case_fnames']
    control_fnames = params['control_fnames']
    sample_fnames = params['sample_fnames']
    gsea_resamples = params['gsea_resamples']
    hg_weight = params['hg_weight']
    min_replace = params['min_replace'],
    noise_mean = params['noise_mean']
    noise_std = params['noise_std']
    pathway_names = params['pathway_names']
    pbar = params['pbar']
    percent = params['percent']
    plage_weight = params['plage_weight']
    prob_missing_peaks = params['prob_missing_peaks']
    random = params['random']

    # constructs the peak intensity dataframe, adding random peaks if necessary
    int_df, updated_pathway_names = construct_intensity_df(sample_fnames, pathway_names, random=random,
                                                           background_pathways=background_pathways)
    int_df = add_random_peaks(sample_fnames, pathway_names, int_df, percent, noise_mean, noise_std)
    ds = convert_to_data_source(int_df, updated_pathway_names, case_fnames, control_fnames, prob_missing_peaks,
                                min_replace)

    # run ORA
    ora = ORA(ds)
    ora_df = ora.get_pathway_df(correct_multiple_tests=True)

    # run PALS
    pals = PALS(ds, plage_weight=plage_weight, hg_weight=hg_weight)
    pals_df = pals.get_pathway_df()

    # run GSEA
    gsea = GSEA(ds, random_sets=gsea_resamples, pbar=pbar)
    gsea_df = gsea.get_pathway_df()

    # store the results
    item = {
        # 'data': ds_resampled, # TOO BIG!!
        'PALS': pals_df,
        'ORA': ora_df,
        'GSEA': gsea_df
    }
    return item


def construct_single_box_df(results, random_peaks, prob_missing, noise_std, method):
    columns = ['pathway', 'percent', 'prob_missing', 'noise_std', 'p_value', 'comb_p_value', 'method']
    box_plot_df = pd.DataFrame(columns=columns)
    try:
        df = pd.concat(results[method], axis=0)
        box_plot_df['pathway'] = list(df.index)
        box_plot_df['percent'] = float(random_peaks)
        box_plot_df['prob_missing'] = float(prob_missing)
        box_plot_df['noise_std'] = float(noise_std)
        box_plot_df['p_value'] = df['case/control p-value'].values
        box_plot_df['comb_p_value'] = df['case/control comb_p'].values
        box_plot_df['method'] = method
        box_plot_df = box_plot_df[~box_plot_df.pathway.str.contains("background")]
        return box_plot_df
    except ValueError:
        return box_plot_df


########################################################################################################################
# Experiments with resampling peaks from real data
########################################################################################################################

def run_resample_experiment(experiment_name, data_source, case, control, proportions, significant_column, n_iter,
                            plage_weight, hg_weight, gsea_resamples, pbar=False, parallel=False):
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
        'hg_weight': hg_weight,
        'gsea_resamples': gsea_resamples
    }

    # vary the number of (mzML) samples, run the pathway analysis methods for n_iter
    n_samples = proportions * data_source.get_measurements().shape[0]
    n_samples = n_samples.astype(int)
    for i in range(len(proportions)):
        # calculate proportions of peaks to resample
        prop = proportions[i]
        n_sample = n_samples[i]
        params = _get_resampled_params(case, control, data_source, gsea_resamples, hg_weight, n_sample,
                                       plage_weight, pbar=pbar)
        logger.info('prop=%.2f n_sample=%d PALS experiment=%s case=%s control=%s' % (
            prop, n_sample, experiment_name, case, control))

        if not parallel:
            prop_results = []  # to store experiment results for this proportion
            for j in range(n_iter):
                item = _get_resampled_performance(params)
                prop_results.append(item)
        else:
            all_params = [params for j in range(n_iter)]
            import ipyparallel as ipp
            rc = ipp.Client()
            dview = rc[:]  # use all engines​
            with dview.sync_imports():
                pass
            prop_results = dview.map_sync(_get_resampled_performance, all_params)

        res['experiment_results'][prop] = prop_results
    return res


def _get_resampled_params(case, control, data_source, gsea_resamples, hg_weight, n_sample, plage_weight, pbar):
    params = {
        'case': case,
        'control': control,
        'data_source': data_source,
        'gsea_resamples': gsea_resamples,
        'hg_weight': hg_weight,
        'n_sample': n_sample,
        'plage_weight': plage_weight,
        'pbar': pbar
    }
    return params


def _get_resampled_performance(params):
    case = params['case']
    control = params['control']
    data_source = params['data_source']
    gsea_resamples = params['gsea_resamples']
    hg_weight = params['hg_weight']
    n_sample = params['n_sample']
    plage_weight = params['plage_weight']
    pbar = params['pbar']

    # generate n_iter data sources, each time randomly drawing n_sample peaks
    ds_resampled = data_source.resample(n_sample, case=case, control=control, axis=0)

    # run ORA on the resampled data
    ora = ORA(ds_resampled, case=case, control=control)
    ora_df = ora.get_pathway_df()

    # run PALS on the resampled data
    pals = PALS(ds_resampled, plage_weight=plage_weight, hg_weight=hg_weight, case=case, control=control)
    pals_df = pals.get_pathway_df()

    # run GSEA on the resampled data
    gsea = GSEA(ds_resampled, random_sets=gsea_resamples, pbar=pbar, case=case, control=control)
    gsea_df = gsea.get_pathway_df()

    # store the results
    item = {
        # 'data': ds_resampled, # TOO BIG!!
        'PALS': pals_df,
        'ORA': ora_df,
        'GSEA': gsea_df
    }
    return item


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
    gsea_resamples = res['gsea_resamples']
    case = res['case']
    control = res['control']
    performances = []

    # generate PALS full results
    pals = PALS(ds, plage_weight=plage_weight, hg_weight=hg_weight, case=case, control=control)
    pals_full_df = pals.get_pathway_df()
    pals_full = _select_significant_entries(pals_full_df, significant_column, SIGNIFICANT_THRESHOLD, N)

    # generate ORA full results
    ora = ORA(ds, case=case, control=control)
    ora_full_df = ora.get_pathway_df()
    ora_full = _select_significant_entries(ora_full_df, significant_column, SIGNIFICANT_THRESHOLD, N)

    # generate GSEA full results
    gsea = GSEA(ds, random_sets=gsea_resamples, pbar=False, case=case, control=control)
    gsea_full_df = gsea.get_pathway_df()
    gsea_full = _select_significant_entries(gsea_full_df, significant_column, SIGNIFICANT_THRESHOLD, N)

    # evaluate the partial results w.r.t to the full results
    logger.debug('Evaluating partial results')
    proportions = list(experiment_results.keys())
    for prop in proportions:
        iterations = range(len(experiment_results[prop]))
        for i in iterations:
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

            # for GSEA
            method = 'GSEA'
            full = gsea_full
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
