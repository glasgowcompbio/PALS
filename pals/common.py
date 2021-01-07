import base64
import copy
import gzip
import json
import logging
import os
import pathlib
import pickle
import random
import sys
from collections import defaultdict
from io import StringIO

import numpy as np
import pandas as pd
from loguru import logger

PW_F_OFFSET = 1
MIN_REPLACE = 5000
MIN_HITS = 0
NUM_RESAMPLES = 1000
PLAGE_WEIGHT = 1
HG_WEIGHT = 0  # TODO: remove this?
SIGNIFICANT_THRESHOLD = 0.05
GSEA_SIGNIFICANT_THRESHOLD = 0.25
SMALL = 1E-6

DATABASE_PIMP_KEGG = 'PiMP_KEGG'
DATABASE_REACTOME_KEGG = 'COMPOUND'
DATABASE_REACTOME_CHEBI = 'ChEBI'
DATABASE_REACTOME_UNIPROT = 'UniProt'
DATABASE_REACTOME_ENSEMBL = 'ENSEMBL'
DATABASE_GNPS_MOLECULAR_FAMILY = 'GNPS'
DATABASE_GNPS_MS2LDA = 'MS2LDA'

GNPS_DOWNLOAD_CYTOSCAPE_DATA_VIEW = 'download_cytoscape_data'
GNPS_VIEW_ALL_MOTIFS_VIEW = 'view_all_motifs'

REACTOME_SPECIES_ARABIDOPSIS_THALIANA = 'Arabidopsis thaliana'
REACTOME_SPECIES_BOS_TAURUS = 'Bos taurus'
REACTOME_SPECIES_CAENORHABDITIS_ELEGANS = 'Caenorhabditis elegans'
REACTOME_SPECIES_CANIS_LUPUS_FAMILIARIS = 'Canis lupus familiaris'
REACTOME_SPECIES_DANIO_RERIO = 'Danio rerio'
REACTOME_SPECIES_DICTYOSTELIUM_DISCOIDEUM = 'Dictyostelium discoideum'
REACTOME_SPECIES_DROSOPHILA_MELANOGASTER = 'Drosophila melanogaster'
REACTOME_SPECIES_GALLUS_GALLUS = 'Gallus gallus'
REACTOME_SPECIES_HOMO_SAPIENS = 'Homo sapiens'
REACTOME_SPECIES_MUS_MUSCULUS = 'Mus musculus'
REACTOME_SPECIES_ORYZA_SATIVA = 'Oryza sativa'
REACTOME_SPECIES_RATTUS_NORVEGICUS = 'Rattus norvegicus'
REACTOME_SPECIES_SACCHAROMYCES_CEREVISIAE = 'Saccharomyces cerevisiae'
REACTOME_SPECIES_SUS_SCROFA = 'Sus scrofa'

PATHWAY_ANALYSIS_PLAGE = 'PLAGE'
PATHWAY_ANALYSIS_ORA = 'ORA'
PATHWAY_ANALYSIS_GSEA = 'GSEA'

GSEA_RANKING_SNR = 'signal_to_noise'
GSEA_RANKING_LOGFC = 'log2_ratio_of_classes'

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
GUI_PATHWAY_ANALYSIS = 'Pathway Analysis'
GUI_GNPS_MOLECULAR_FAMILY_ANALYSIS = 'GNPS Molecular Family Analysis'
GUI_GNPS_MS2LDA_ANALYSIS = 'GNPS-MS2LDA Analysis'


def load_json(json_file, compressed=False):
    if compressed:
        with gzip.GzipFile(json_file, 'r') as f:
            data = json.loads(f.read().decode('utf-8'))
    else:
        with open(json_file, 'r') as f:
            data = json.load(f)
    return data


def save_json(data, json_file, compressed=False):
    if compressed:
        with gzip.GzipFile(json_file, 'w') as f:
            f.write(json.dumps(data).encode('utf-8'))
    else:
        with open(json_file, 'w') as f:
            json.dump(data, f)


def create_if_not_exist(out_dir):
    if not os.path.exists(out_dir) and len(out_dir) > 0:
        logger.debug('Created %s' % out_dir)
        pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)


def save_obj(obj, filename):
    """
    Save object to file
    :param obj: the object to save
    :param filename: the output file
    :return: None
    """
    out_dir = os.path.dirname(filename)
    create_if_not_exist(out_dir)
    logger.debug('Saving %s to %s' % (type(obj), filename))
    with gzip.GzipFile(filename, 'w') as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_obj(filename):
    """
    Load saved object from file
    :param filename: The file to load
    :return: the loaded object
    """
    try:
        with gzip.GzipFile(filename, 'rb') as f:
            return pickle.load(f)
    except OSError:
        logger.warning('Old, invalid or missing pickle in %s. Being regenerated.' % filename)
        return None


def set_log_level_warning():
    logger.remove()
    logger.add(sys.stderr, level=logging.WARNING)


def set_log_level_info():
    logger.remove()
    logger.add(sys.stderr, level=logging.INFO)


def set_log_level_debug():
    logger.remove()
    logger.add(sys.stderr, level=logging.DEBUG)


def is_comparison_used(comp, selected_case, selected_control):
    # skip this comparison if not the same as what the user has specified
    comp_case = comp['case']
    comp_control = comp['control']
    if selected_case is not None and selected_control is not None:
        if selected_case != comp_case or selected_control != comp_control:
            return False
    return True


def load_data(intensity_csv, annotation_csv, gui=False):
    """
    Loads PALS data from csv files
    :param intensity_csv: a CSV file of peak intensities
    :param annotation_csv: a CSV file of peak annotation
    :return: the peak intensity df, annotation df, and also grouping information
    """
    # load intensity dataframe from csv
    # skip the second row containing group information
    # set the first column (peak ids) in csv to be the index

    if gui:  # load data for GUI (Streamlit)
        lines = uploaded_file_to_strings(intensity_csv)  # turn to list of strings
        filtered_lines = [line for idx, line in enumerate(lines) if idx != 1]  # remove second (grouping) line
        intensities = StringIO('\n'.join(filtered_lines))

        # load to pandas dataframe (without the second line)
        int_df = pd.read_csv(intensities, skiprows=[1], index_col=0)
        groups = _parse_groups(lines)
        logger.debug('Loaded %d x %d peak intensities from %s' % (int_df.shape[0], int_df.shape[1], intensity_csv))
        logger.debug('Loaded groups: %s' % groups)

        # load annotation dataframe from csv, setting the first column to be the index
        lines = uploaded_file_to_strings(annotation_csv)
        annotations = StringIO('\n'.join(lines))
        annotation_df = pd.read_csv(annotations, index_col=0)
        logger.debug('Loaded %d peak annotations from %s' % (annotation_df.shape[0], annotation_csv))

    else:  # for command-line use
        int_df = pd.read_csv(intensity_csv, skiprows=[1], index_col=0)
        groups = get_groups(intensity_csv)
        logger.debug('Loaded %d x %d peak intensities from %s' % (int_df.shape[0], int_df.shape[1], intensity_csv))
        logger.debug('Loaded groups: %s' % groups)

        annotation_df = pd.read_csv(annotation_csv, index_col=0)
        logger.debug('Loaded %d peak annotations from %s' % (annotation_df.shape[0], annotation_csv))

    return int_df, annotation_df, groups


def uploaded_file_to_strings(csv_file):
    bytesData = csv_file.getvalue()
    stringio = StringIO(str(bytesData, 'utf-8'))
    lines = stringio.readlines()
    return lines


def get_groups(csv_file):
    """
    Extract group information from peak intensity csv file
    :param csv_file: An input CSV file.
    The first row in the CSV file should contain samples information.
    The second row in the CSV file should contain grouping information.
    :return: a dictionary where key is the group name and values are the samples under that group
    """
    with open(csv_file, 'r') as f:
        lines = f.readlines()
        grouping = _parse_groups(lines)
        return grouping


def _parse_groups(lines):
    grouping = defaultdict(list)
    # extract the first and second lines containing samples and grouping information
    samples = None
    groups = None
    for i, line in enumerate(lines):
        if i == 0:
            samples = line.strip().split(',')
        if i == 1:
            groups = line.strip().split(',')
            break
    assert samples is not None, 'Missing samples line'
    assert groups is not None, 'Missing groups line'
    for sample, group in zip(samples, groups):
        if group == 'group':  # skip the column header containing the word 'group'
            continue
        grouping[group].append(sample)
    return dict(grouping)


class Method(object):
    def __init__(self, data_source, seed=None, preprocessors=None):
        if seed is None:
            random.seed()
        else:
            random.seed(seed)

        self.data_source = copy.deepcopy(data_source)
        if preprocessors is None:
            self.preprocessors = self._create_preprocessors()

    def get_results(self, preprocess=True):
        raise NotImplementedError()

    def get_pathway_df(self, standardize=True): # to support old method calls
        return self.get_results(preprocess=standardize)

    def _create_preprocessors(self):
        return []

    def _get_measurement_df(self, preprocess):
        df = self.data_source.get_measurements()
        if preprocess:
            for p in self.preprocessors:
                df = p.process(df)
        return df


# https://discuss.streamlit.io/t/how-to-set-file-download-function/2141
def get_table_download_link(df, out_file, label):
    """Generates a link allowing the data in a given panda dataframe to be downloaded
    in:  dataframe
    out: href string
    """
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(
        csv.encode()
    ).decode()  # some strings <-> bytes conversions necessary here
    out = f'<a href="data:file/csv;base64,{b64}" download="%s">%s</a>' % (out_file, label)
    return out


def post_filter_df_by_min_hits(pathway_df, min_hits):
    """
    Filter the 'tot_ds_F' column of a pathway dataframe. If it's less than min_hits, then set all the p-values to NaN.
    :param pathway_df: a pathway dataframe returned by a pathway ranking method in PALS
    :param min_hits: the minimum number of hits in this pathway
    :return: a filtered pathway dataframe
    """
    pathway_df = pathway_df.copy()
    rows_to_remove = pathway_df[pathway_df['tot_ds_F'] < min_hits]
    columns_to_remove = [col for col in rows_to_remove.columns.values if 'p-value' in col or 'comb_p' in col]
    pathway_df.loc[rows_to_remove.index, columns_to_remove] = np.nan
    return pathway_df
