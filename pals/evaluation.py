import warnings

import numpy as np
import pandas as pd
from loguru import logger
from matplotlib.patches import PathPatch
from sklearn.metrics import auc

from .GSEA import GSEA
from .ORA import ORA
from .PLAGE import PLAGE
from .common import SIGNIFICANT_THRESHOLD, MIN_REPLACE, set_log_level_info, set_log_level_debug, NUM_RESAMPLES, \
    PLAGE_WEIGHT, HG_WEIGHT, GSEA_RANKING_SNR, GSEA_SIGNIFICANT_THRESHOLD
from .noise import construct_intensity_df, add_random_peaks, convert_to_data_source


########################################################################################################################
# Experiments with adding noise on synthetic data
########################################################################################################################

def run_noise_experiment(background_pathways, case_fnames, control_fnames, pathway_names, num_iterations, plage_weight,
                         hg_weight, gsea_resamples, gsea_ranking_method, reqd_scenarios, parallel=False):
    pals_dfs = []
    ora_dfs = []
    gsea_dfs = []
    random = False  # whether to sample intensity data randomly or using pre-set values

    exp_results = {}
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
                                       hg_weight=hg_weight, gsea_resamples=gsea_resamples,
                                       gsea_ranking_method=gsea_ranking_method, parallel=parallel)

        df = construct_single_box_df(results, scenario['percent'], scenario['prob_missing_peaks'],
                                     scenario['noise_std'], 'PALS')
        pals_dfs.append(df)

        df = construct_single_box_df(results, scenario['percent'], scenario['prob_missing_peaks'],
                                     scenario['noise_std'], 'ORA')
        ora_dfs.append(df)

        df = construct_single_box_df(results, scenario['percent'], scenario['prob_missing_peaks'],
                                     scenario['noise_std'], 'GSEA')
        gsea_dfs.append(df)
        exp_results[i] = results

    pals_df = pd.concat(pals_dfs, axis=0)
    ora_df = pd.concat(ora_dfs, axis=0)
    gsea_df = pd.concat(gsea_dfs, axis=0)

    return pals_df, ora_df, gsea_df, exp_results


def calc_av_p_scores(case_fnames, control_fnames, pathway_names, num_iterations, percent=0, random=False, noise_mean=0,
                     noise_std=5, prob_missing_peaks=0.2, background_pathways=100, min_replace=MIN_REPLACE,
                     plage_weight=PLAGE_WEIGHT, hg_weight=HG_WEIGHT, gsea_resamples=NUM_RESAMPLES,
                     gsea_ranking_method=GSEA_RANKING_SNR, parallel=False):
    sample_fnames = control_fnames + case_fnames
    set_log_level_info()

    # create num_iterations random data
    params = _get_noise_params(background_pathways, case_fnames, control_fnames, gsea_resamples, gsea_ranking_method,
                               hg_weight, min_replace, noise_mean, noise_std, pathway_names, percent, plage_weight,
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


def _get_noise_params(background_pathways, case_fnames, control_fnames, gsea_resamples, gsea_ranking_method, hg_weight,
                      min_replace, noise_mean, noise_std, pathway_names, percent, plage_weight,
                      prob_missing_peaks, random, sample_fnames):
    return {
        'background_pathways': background_pathways,
        'case_fnames': case_fnames,
        'control_fnames': control_fnames,
        'gsea_resamples': gsea_resamples,
        'gsea_ranking_method': gsea_ranking_method,
        'hg_weight': hg_weight,
        'min_replace': min_replace,
        'noise_mean': noise_mean,
        'noise_std': noise_std,
        'pathway_names': pathway_names,
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
    gsea_ranking_method = params['gsea_ranking_method']
    hg_weight = params['hg_weight']
    min_replace = params['min_replace'],
    noise_mean = params['noise_mean']
    noise_std = params['noise_std']
    pathway_names = params['pathway_names']
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
    pals = PLAGE(ds, plage_weight=plage_weight, hg_weight=hg_weight)
    pals_df = pals.get_pathway_df()

    # run GSEA
    gsea = GSEA(ds, num_resamples=gsea_resamples, method=gsea_ranking_method)
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
        return box_plot_df
    except ValueError:
        return box_plot_df


def get_tp_fn_fn(reqd_scenarios, exp_results, true_answers):
    data = []
    for i in exp_results:
        scenario = reqd_scenarios[i]
        noise_std = scenario['noise_std']
        percent = scenario['percent']
        prob_missing_peaks = scenario['prob_missing_peaks']

        results = exp_results[i]
        for method in results:
            dataframes = results[method]
            for df in dataframes:
                if method == 'GSEA':
                    threshold = SIGNIFICANT_THRESHOLD
                else:
                    threshold = SIGNIFICANT_THRESHOLD
                filtered_df = _select_significant_entries(df, 'case/control comb_p', threshold, None)
                TP, FP, FN, prec, rec, f1, TP_items, FP_items, FN_items = _compute_prec_rec_f1(
                    true_answers, set(filtered_df.index.values))
                row = [method, noise_std, percent, prob_missing_peaks, TP, FP, FN, prec, rec, f1, TP_items, FP_items,
                       FN_items]
                data.append(row)

    df = pd.DataFrame(data,
                      columns=['method', 'noise_std', 'percent', 'prob_missing_peaks', 'TP', 'FP', 'FN', 'prec', 'rec',
                               'f1', 'TP_items', 'FP_items', 'FN_items'])
    return df


def compute_pr_curve(method, df, true_answers):
    pr_results = []
    for threshold in df['p_value'].unique():
        filtered_df = df[df['p_value'] <= threshold]
        TP, FP, FN, prec, rec, f1, TP_items, FP_items, FN_items = _compute_prec_rec_f1(true_answers,
                                                                                       set(filtered_df.index.values))
        row = (method, threshold, TP, FP, FN, prec, rec, f1)
        pr_results.append(row)

    pr_df = pd.DataFrame(pr_results, columns=['method', 'threshold', 'TP', 'FP', 'FN', 'prec', 'rec', 'f1'])
    sorted_pr_df = pr_df.sort_values('prec', ascending=False)
    return sorted_pr_df


def get_auc_for_noise(reqd_scenarios, exp_results, true_answers):
    significant_column = 'case/control comb_p'
    auc_results = []
    for i in range(len(reqd_scenarios)):
        scenario = reqd_scenarios[i]
        results = exp_results[i]

        for method in results:
            dataframes = results[method]
            method_aucs = []

            for j in range(len(dataframes)):
                pathway_df = dataframes[j]
                df = pathway_df.sort_values(significant_column, ascending=False)
                df = df.rename(columns={significant_column: 'p_value'})
                try:
                    df = df[['pw_name', 'p_value', 'sf', 'unq_pw_F', 'tot_ds_F', 'F_coverage']]
                except KeyError:
                    df = df[['pw_name', 'p_value', 'unq_pw_F', 'tot_ds_F', 'F_coverage']]

                try:
                    sorted_pr_df = compute_pr_curve(method, df, true_answers)
                    auc_res = auc(sorted_pr_df['prec'], sorted_pr_df['rec'])  # compute auc
                except ValueError:
                    auc_res = 0.0
                method_aucs.append(auc_res)

            # show progress by printing average auc
            avg_method_auc = np.mean(method_aucs)
            print('%s %s %.3f' % (scenario, method, avg_method_auc))

            for j in range(len(method_aucs)):
                auc_res = method_aucs[j]
                auc_results.append(
                    (method, scenario['noise_std'], scenario['percent'], scenario['prob_missing_peaks'], j, auc_res))

        print()

    auc_df = pd.DataFrame(auc_results, columns=['method', 'noise_std', 'percent', 'prob_missing_peaks', 'iter', 'auc'])
    return auc_df


########################################################################################################################
# Experiments with resampling peaks from real data
########################################################################################################################

def run_resample_experiment(experiment_name, data_source, case, control, prob_missing_peaks, significant_column, n_iter,
                            plage_weight, hg_weight, gsea_resamples, gsea_ranking_method, parallel=False):
    res = {
        'experiment_name': experiment_name,
        'data_source': data_source,
        'case': case,
        'control': control,
        'prob_missing_peaks': prob_missing_peaks,
        'n_iter': n_iter,
        'significant_column': significant_column,
        'experiment_results': {},
        'plage_weight': plage_weight,
        'hg_weight': hg_weight,
        'gsea_resamples': gsea_resamples,
        'gsea_ranking_method': gsea_ranking_method
    }

    # vary the number of (mzML) samples, run the pathway analysis methods for n_iter
    n_samples = (1-prob_missing_peaks) * data_source.get_measurements().shape[0]
    n_samples = n_samples.astype(int)
    for i in range(len(prob_missing_peaks)):
        # calculate proportions of peaks to resample
        prob = prob_missing_peaks[i]
        n_sample = n_samples[i]
        params = _get_resampled_params(case, control, data_source, gsea_resamples, gsea_ranking_method, hg_weight,
                                       n_sample, plage_weight)
        logger.info('prob_missing_peaks=%.2f n_sample=%d PALS experiment=%s case=%s control=%s' % (
            prob, n_sample, experiment_name, case, control))

        if not parallel:
            prob_results = []  # to store experiment results for this prob value
            for j in range(n_iter):
                item = _get_resampled_performance(params)
                prob_results.append(item)
        else:
            all_params = [params for j in range(n_iter)]
            import ipyparallel as ipp
            rc = ipp.Client()
            dview = rc[:]  # use all engines​
            with dview.sync_imports():
                pass
            prob_results = dview.map_sync(_get_resampled_performance, all_params)

        res['experiment_results'][prob] = prob_results
    return res


def _get_resampled_params(case, control, data_source, gsea_resamples, gsea_ranking_method, hg_weight, n_sample,
                          plage_weight):
    params = {
        'case': case,
        'control': control,
        'data_source': data_source,
        'gsea_resamples': gsea_resamples,
        'gsea_ranking_method': gsea_ranking_method,
        'hg_weight': hg_weight,
        'n_sample': n_sample,
        'plage_weight': plage_weight,
    }
    return params


def _get_resampled_performance(params):
    case = params['case']
    control = params['control']
    data_source = params['data_source']
    gsea_resamples = params['gsea_resamples']
    gsea_ranking_method = params['gsea_ranking_method']
    hg_weight = params['hg_weight']
    n_sample = params['n_sample']
    plage_weight = params['plage_weight']

    # generate n_iter data sources, each time randomly drawing n_sample peaks
    ds_resampled = data_source.resample(n_sample, case=case, control=control, axis=0)

    # run ORA on the resampled data
    ora = ORA(ds_resampled, case=case, control=control)
    ora_df = ora.get_pathway_df()

    # run PALS on the resampled data
    pals = PLAGE(ds_resampled, plage_weight=plage_weight, hg_weight=hg_weight, case=case, control=control)
    pals_df = pals.get_pathway_df()

    # run GSEA on the resampled data
    gsea = GSEA(ds_resampled, num_resamples=gsea_resamples, method=gsea_ranking_method, case=case, control=control)
    gsea_df = gsea.get_pathway_df()

    # store the results
    item = {
        # 'data': ds_resampled, # TOO BIG!!
        'PALS': pals_df,
        'ORA': ora_df,
        'GSEA': gsea_df
    }
    return item


def get_method_true_answers(res, N=None):
    significant_column = res['significant_column']
    ds = res['data_source']
    plage_weight = res['plage_weight']
    hg_weight = res['hg_weight']
    gsea_resamples = res['gsea_resamples']
    gsea_ranking_method = res['gsea_ranking_method']
    case = res['case']
    control = res['control']

    # generate the true answers for each method
    pals = PLAGE(ds, plage_weight=plage_weight, hg_weight=hg_weight, case=case, control=control)
    pals_full_df = pals.get_pathway_df()
    pals_full = _select_significant_entries(pals_full_df, significant_column, SIGNIFICANT_THRESHOLD, N)

    ora = ORA(ds, case=case, control=control)
    ora_full_df = ora.get_pathway_df()
    ora_full = _select_significant_entries(ora_full_df, significant_column, SIGNIFICANT_THRESHOLD, N)

    gsea = GSEA(ds, num_resamples=gsea_resamples, method=gsea_ranking_method, case=case, control=control)
    gsea_full_df = gsea.get_pathway_df()
    gsea_full = _select_significant_entries(gsea_full_df, significant_column, GSEA_SIGNIFICANT_THRESHOLD, N)

    method_true_answers = {
        'PALS': pals_full,
        'ORA': ora_full,
        'GSEA': gsea_full,
    }
    return method_true_answers


def evaluate_performance(res, true_answers, N=None):
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
    experiment_results = res['experiment_results']
    significant_column = res['significant_column']
    pals_full = true_answers['PALS']
    ora_full = true_answers['ORA']
    gsea_full = true_answers['GSEA']

    # evaluate the partial results w.r.t to the full results
    performances = []
    logger.debug('Evaluating partial results')
    prob_missing_peaks = list(experiment_results.keys())
    for prob in prob_missing_peaks:
        iterations = range(len(experiment_results[prob]))
        for i in iterations:
            item = experiment_results[prob][i]

            # for PALS
            method = 'PALS'
            full = pals_full
            partial_df = item[method]
            partial = _select_significant_entries(partial_df, significant_column, SIGNIFICANT_THRESHOLD, N)
            TP, FP, FN, prec, rec, f1, TP_items, FP_items, FN_items = _compute_prec_rec_f1(set(full.index.values),
                                                                                           set(partial.index.values))
            perf = (method, prob, i, TP, FP, FN, prec, rec, f1)
            performances.append(perf)

            # for ORA
            method = 'ORA'
            full = ora_full
            partial_df = item[method]
            partial = _select_significant_entries(partial_df, significant_column, SIGNIFICANT_THRESHOLD, N)
            TP, FP, FN, prec, rec, f1, TP_items, FP_items, FN_items = _compute_prec_rec_f1(set(full.index.values),
                                                                                           set(partial.index.values))
            perf = (method, prob, i, TP, FP, FN, prec, rec, f1)
            performances.append(perf)

            # for GSEA
            method = 'GSEA'
            full = gsea_full
            partial_df = item[method]
            partial = _select_significant_entries(partial_df, significant_column, GSEA_SIGNIFICANT_THRESHOLD, N)
            TP, FP, FN, prec, rec, f1, TP_items, FP_items, FN_items = _compute_prec_rec_f1(set(full.index.values),
                                                                                           set(partial.index.values))
            perf = (method, prob, i, TP, FP, FN, prec, rec, f1)
            performances.append(perf)

    logger.debug('Done!')
    performance_df = pd.DataFrame(performances,
                                  columns=['method', 'missing_peaks', 'i', 'TP', 'FP', 'FN', 'precision', 'recall', 'F1'])
    return performance_df


def get_auc_for_hat_data(res, true_answers):
    aucs = []
    exp_res = res['experiment_results']
    significant_column = res['significant_column']

    for prob in exp_res:
        logger.debug('Processing %f' % prob)
        prob_res = exp_res[prob]
        for i in range(len(prob_res)):
            method_res = prob_res[i]
            for method in method_res:
                pathway_df = method_res[method]

                df = pathway_df.sort_values(significant_column, ascending=False)
                df = df.rename(columns={significant_column: 'p_value'})
                try:
                    df = df[['pw_name', 'p_value', 'sf', 'unq_pw_F', 'tot_ds_F', 'F_coverage']]
                except KeyError:
                    df = df[['pw_name', 'p_value', 'unq_pw_F', 'tot_ds_F', 'F_coverage']]

                true_answer_df = true_answers[method]
                true_answer_set = set(true_answer_df.index.values.tolist())
                sorted_pr_df = compute_pr_curve(method, df, true_answer_set)
                computed_auc = auc(sorted_pr_df['prec'], sorted_pr_df['rec'])  # compute auc
                aucs.append((prob, i, method, computed_auc))

    auc_df = pd.DataFrame(aucs, columns=['missing_peaks', 'i', 'method', 'auc'])
    return auc_df


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


def _compute_prec_rec_f1(pathways_full, pathways_partial):
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
    return TP, FP, FN, prec, rec, f1, TP_items, FP_items, FN_items


# https://stackoverflow.com/questions/56838187/how-to-create-spacing-between-same-subgroup-in-seaborn-boxplot
def adjust_box_widths(g, fac):
    """
    Adjust the withs of a seaborn-generated boxplot.
    """

    # iterating through Axes instances
    for ax in g.axes:

        # iterating through axes artists:
        for c in ax.get_children():

            # searching for PathPatches
            if isinstance(c, PathPatch):
                # getting current width of box:
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5 * (xmin + xmax)
                xhalf = 0.5 * (xmax - xmin)

                # setting new width of box
                xmin_new = xmid - fac * xhalf
                xmax_new = xmid + fac * xhalf
                verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new

                # setting new width of median line
                for l in ax.lines:
                    if np.all(l.get_xdata() == [xmin, xmax]):
                        l.set_xdata([xmin_new, xmax_new])
