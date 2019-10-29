import os
import pickle

import pandas as pd

from feature_extraction import DataSource
from pathway_analysis import PALS
from pimp_tools import get_pimp_API_token_from_env, PIMP_HOST, get_ms1_intensities, get_formula_df, \
    get_experimental_design

token = get_pimp_API_token_from_env()
analysis_id = 1321  # example beer analysis

int_df_filename = os.path.join(os.getcwd(), 'notebooks', 'test_data', 'int_df.p')
try:
    int_df = pd.read_pickle(int_df_filename)
except FileNotFoundError:
    int_df = get_ms1_intensities(token, PIMP_HOST, analysis_id)
    int_df.to_pickle(int_df_filename)

formula_df_filename = os.path.join(os.getcwd(), 'notebooks', 'test_data', 'formula_df.p')
try:
    formula_df = pd.read_pickle(formula_df_filename)
except FileNotFoundError:
    formula_df = get_formula_df(token, PIMP_HOST, analysis_id)
    formula_df.to_pickle('formula_df.p')

experimental_design_filename = os.path.join(os.getcwd(), 'test_data', 'experimental_design.p')
try:
    experimental_design_filename = os.path.join(os.getcwd(), 'test_data', 'experimental_design.p')
    with open(experimental_design_filename, 'rb') as f:
        experimental_design = pickle.load(f)
except FileNotFoundError:
    experimental_design = get_experimental_design(token, PIMP_HOST, analysis_id)
    with open(experimental_design_filename, 'wb') as f:
        pickle.dump(experimental_design, f)

ds = DataSource(int_df, formula_df, experimental_design, database_name='kegg')
pals = PALS(ds, min_intensity=5000, num_resamples=10)
activity_df = pals.get_plage_activity_df()
plage_df = pals.set_up_resample_plage_p_df(activity_df)
pathway_df = pals.calculate_hg_values(plage_df)
