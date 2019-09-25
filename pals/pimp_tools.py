import json
import os

import pandas as pd
import requests

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
        print(url, r)

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
    return payload


def get_ms2_peaks(token, host, analysis_id, as_dataframe=False):
    url = 'http://{}/export/get_ms2_peaks?analysis_id={}&as_dataframe={}'.format(host, analysis_id, as_dataframe)
    payload = get_data(token, url, as_dataframe)
    return payload


def get_formula_df(token, host, analysis_id):
    ms1_df = get_ms1_peaks(token, host, analysis_id)
    ms1_df['identified'] = ms1_df['identified'].astype('bool') # convert identified column ('True', 'False') to boolean

    # select only peaks that have been (identified) or (annotated with adduct type M+H and M-H).
    identified_peaks = ms1_df.query('identified == True')
    annotated_peaks = ms1_df.query('identified == False & (adduct == "M+H" | adduct == "M-H")')
    identified_annotated_peaks = pd.concat([identified_peaks, annotated_peaks])

    # set pid as index, extract formula
    formula_df = identified_annotated_peaks[['pid', 'db', 'identifier', 'formula']]
    formula_df = formula_df.set_index('pid')
    return formula_df
