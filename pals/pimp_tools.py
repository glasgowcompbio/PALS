import json
import os
import tempfile

import numpy as np
import pandas as pd
import requests
from loguru import logger

from pals.common import load_obj, save_obj

PIMP_HOST = 'polyomics.mvls.gla.ac.uk'


def get_pimp_API_token_from_env():
    return os.environ['PIMP_API_TOKEN']


def get_authentication_token(host, username, password):
    url = 'http://{}/export/get_token'.format(host)
    r = requests.post(url, data={'username': username, 'password': password})
    token = json.loads(r.text)['token']
    return (token)


def get_data(token, url, as_dataframe=False):
    headers = {'Authorization': 'token {}'.format(token)}
    payload = None

    with requests.Session() as s:
        s.headers.update(headers)
        r = s.get(url)
        logger.debug(url, r)

        # extract GET results
        if r.status_code == 200:
            payload = json.loads(r.text)
            if as_dataframe:
                try:
                    df = pd.read_json(payload)
                except:  # alternative way to load the response
                    df = pd.read_json(r.text)
                payload = df.sort_index()

    return payload


def get_ms1_peaks(token, host, analysis_id):
    url = 'http://{}/export/get_ms1_peaks?analysis_id={}'.format(host, analysis_id)
    payload = get_data(token, url, True)
    return payload


def get_ms1_intensities(token, host, analysis_id):
    url = 'http://{}/export/get_ms1_intensities?analysis_id={}'.format(host, analysis_id)
    payload = get_data(token, url, True)
    payload.index.name = 'row_id'
    return payload


def get_ms2_peaks(token, host, analysis_id, as_dataframe=False):
    url = 'http://{}/export/get_ms2_peaks?analysis_id={}&as_dataframe={}'.format(host, analysis_id, as_dataframe)
    payload = get_data(token, url, as_dataframe)
    return payload


def get_experimental_design(token, host, analysis_id):
    url = 'http://{}/export/get_experimental_design?analysis_id={}'.format(host, analysis_id)
    payload = get_data(token, url, False)
    return payload


def get_annotation_df(token, host, analysis_id, database_name='kegg', polarity='positive'):
    ms1_df = get_ms1_peaks(token, host, analysis_id)
    ms1_df['identified'] = ms1_df['identified'].astype('bool')  # convert identified column ('True', 'False') to boolean
    ms1_df = ms1_df[ms1_df['db'] == database_name]  # filter by db name, e.g. 'kegg'
    # ms1_df = ms1_df[ms1_df['polarity'] == polarity] # filter by polarity, e.g. 'positive'

    ms1_df.rename(columns={'pid': 'row_id', 'identifier': 'entity_id'}, inplace=True)
    ms1_df = ms1_df.set_index('row_id')

    # select only peaks that have been (identified) or (annotated with adduct type M+H and M-H).
    # doesn't work
    # identified_peaks = ms1_df.query('identified == True')
    # annotated_peaks = ms1_df.query('identified == False & (adduct == "M+H" | adduct == "M-H")')
    # identified_annotated_peaks = pd.concat([identified_peaks, annotated_peaks])
    # formula_df = identified_annotated_peaks

    # from PiMP export, 'identified' is somehow always true, i.e. above logic doesn't work
    # below is an alternative implementation of the same filtering
    formulas = ms1_df['formula'].unique()
    to_remove = []
    for formula in formulas:
        peaks = ms1_df[ms1_df['formula'] == formula]

        # filter by adducts
        adducts = peaks['adduct'].unique()
        no_adduct = 'M+H' not in adducts and 'M-H' not in adducts

        # filter by ms2 annotations
        frank_annots = peaks['frank_annot'].values
        no_annot = np.all(pd.isnull(frank_annots))

        # if peaks not M+H or M-H and don't have ms2 annotation, then add the formula for removal
        if no_adduct and no_annot:
            to_remove.append(formula)

    # keep peaks having formulas NOT in the to_remove list
    formula_df = ms1_df[~ms1_df['formula'].isin(to_remove)]

    # select only the column we need
    formula_df = formula_df[['entity_id']]
    return formula_df


def download_from_pimp(token, hostname, analysis_id, database_name):
    temp_dir = tempfile.gettempdir()
    filename = os.path.join(temp_dir, 'pimp_analysis_%d.p' % analysis_id)
    logger.debug('Trying to load data from temp file: %s' % filename)

    data = load_obj(filename)
    if data is None:
        logger.debug('Retrieving data for analysis %d from PiMP' % analysis_id)
        int_df = get_ms1_intensities(token, hostname, analysis_id)
        annotation_df = get_annotation_df(token, hostname, analysis_id, database_name)
        experimental_design = get_experimental_design(token, hostname, analysis_id)
        data = {
            'int_df': int_df,
            'annotation_df': annotation_df,
            'experimental_design': experimental_design
        }
        logger.debug('Caching analysis data for next use')
        save_obj(data, filename)
    else:
        int_df = data['int_df']
        annotation_df = data['annotation_df']
        experimental_design = data['experimental_design']
    return int_df, annotation_df, experimental_design
