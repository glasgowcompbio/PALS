import math
from collections import defaultdict

import numpy as np
import pandas as pd
import pylab as plt
import seaborn as sns
from loguru import logger

from pals.feature_extraction import DataSource


def construct_intensity_df(sample_fnames, pathway_names, random=False, background_pathways=None):
    # some predetermined data with fold changes between conditions
    synthetic_data = [20.0, 20.0, 20.0, 20.0, 40.0, 40.0, 40.0, 40.0]

    # randomly sample some noise from a normal distribution
    random_data = np.random.normal(0, 1, len(sample_fnames))

    # copy pathway names dictionary
    pathway_names = dict(pathway_names)

    # loop over pathway names, generate logged peak data (with some noise) up to num peaks
    pk_samp_intensities = []
    for name, num in pathway_names.items():
        for n in range(num):
            peak_int_list = []
            peak_int_list.append(name)
            if not random:
                data_noise = synthetic_data + np.random.normal(0, 5, len(synthetic_data))
            else:
                data_noise = random_data + np.random.normal(0, 5, len(random_data))
            peak_int_list.extend(list(data_noise))  # The intensities of all the samples for this peak.
            pk_samp_intensities.append(peak_int_list)

    # generate background pathways randomly
    if background_pathways is not None:
        for i in range(background_pathways):
            # sample number of entities in the pathway randomly from 5 .. 50
            num = np.random.randint(5, 51)
            name = 'background%d' % (i)
            pathway_names[name] = num
            for n in range(num):
                peak_int_list = []
                peak_int_list.append(name)
                data_noise = random_data + np.random.normal(0, 5, len(random_data))
                peak_int_list.extend(list(data_noise))  # The intensities of all the samples for this peak.
                pk_samp_intensities.append(peak_int_list)

    int_df = pd.DataFrame(pk_samp_intensities).set_index([0])
    int_df.columns = sample_fnames
    int_df.index.name = "ms1_peak_id"
    int_df.columns.name = "sample_name"
    return int_df, pathway_names


def add_random_peaks(sample_fnames, pathway_names, int_df, percent, noise_mean, noise_std):
    if percent == 0:
        new_df = int_df.copy()
        return new_df
    else:
        # For each of the pathways add random peaks
        rand_peak_list = []
        for name in pathway_names:
            # get all peak intensities in that pathway
            df_path = int_df.loc[name]

            # compute percentage of random peaks to add
            num_peaks = math.ceil((df_path.shape[0]) * (percent / 100.0))
            num_samples = df_path.shape[1]

            # generate random peaks for that pathway
            for p in range(int(num_peaks)):
                rand_peaks = []
                data = np.random.normal(noise_mean, noise_std, num_samples)
                rand_peaks.append(name)
                rand_peaks.extend(list(data))
                rand_peak_list.append(rand_peaks)

        # construct a DF for the new peak list
        ran_df = pd.DataFrame(rand_peak_list).set_index([0])
        ran_df.columns = sample_fnames

        # add the random peaks to the original DF
        new_df = pd.concat([int_df, ran_df])
        new_df.index.name = int_df.index.name
        return new_df


def plot_intensity_matrix(int_df, out_file=None):
    plt.figure(figsize=(5, 10))
    sns.heatmap(int_df, cmap="rainbow")
    plt.title("Simulated intensity matrix (log)")
    plt.tight_layout()
    if out_file is not None:
        plt.savefig(out_file, dpi=300)


def convert_to_data_source(int_df, pathway_names, case_fnames, control_fnames, prob_missing_peaks, min_replace):
    int_df = np.exp(int_df.copy())
    int_df = int_df.reset_index()
    int_df.index = np.arange(1, len(int_df) + 1)
    int_df.index.name = 'row_id'
    int_df = int_df.rename(columns={'ms1_peak_id': 'pathway_id'})
    data, annotations = _get_database(int_df, pathway_names)

    int_df = int_df.drop('pathway_id', axis=1)
    annotation_df = pd.DataFrame.from_dict(annotations, orient='index')
    annotation_df.index.name = 'row_id'
    annotation_df.columns = ['entity_id']
    logger.debug('Dataset annotations = %d' % annotation_df.shape[0])

    # randomly sample (1-prob_missing_peaks) rows from annotation_df without replacement
    # annotation_df = annotation_df.sample(frac=(1 - prob_missing_peaks), replace=False)
    logger.debug('Sampled annotations = %d with prob_missing_peaks=%.2f' % (annotation_df.shape[0], prob_missing_peaks))

    experimental_design = {
        'comparisons': [
            {'case': 'case', 'control': 'control', 'name': 'case/control'}
        ],
        'groups': {
            'case': case_fnames,
            'control': control_fnames
        }
    }
    ds = DataSource(int_df, annotation_df, experimental_design, None, database=data, min_replace=min_replace)
    n_sample = int((1 - prob_missing_peaks) * int_df.shape[0])
    ds_resampled = ds.resample(n_sample, axis=0)
    return ds_resampled


def _get_database(int_df, pathway_names):
    # create pathway dict
    pathway_dict = {}
    for k, v in pathway_names.items():
        pathway_dict[k] = {'display_name': k}

    # create entity and mapping dict
    entity_dict = {}
    mapping_dict = defaultdict(list)
    annotations = {}
    for row_id, row in int_df.iterrows():
        compound_id = 'C%d' % row_id
        pathway_id = row['pathway_id']
        entity_dict[compound_id] = {'unique_id': compound_id, 'display_name': compound_id}
        mapping_dict[compound_id].append(pathway_id)
        annotations[row_id] = compound_id
    mapping_dict = dict(mapping_dict)

    data = {
        'pathway_dict': pathway_dict,
        'entity_dict': entity_dict,
        'mapping_dict': mapping_dict
    }
    return data, annotations
