import pickle
import pandas as pd

# %%

from pals.pimp_tools import get_pimp_API_token_from_env, PIMP_HOST, get_ms1_intensities, get_ms1_peaks, get_formula_df, \
    get_experimental_design
from pals.feature_extraction import DataSource
from pals.pathway_analysis import PALS
from pals.common import *

token = get_pimp_API_token_from_env()
analysis_id = 1321  # example beer analysis

int_df_filename = os.path.join(os.getcwd(), 'test_data', 'int_df.p')
try:
    int_df = pd.read_pickle(int_df_filename)
except FileNotFoundError:
    int_df = get_ms1_intensities(token, PIMP_HOST, analysis_id)
    int_df.to_pickle(int_df_filename)

database_name = 'kegg'
formula_df_filename = os.path.join(os.getcwd(), 'test_data', 'formula_df.p')
try:
    formula_df = pd.read_pickle(formula_df_filename)
except FileNotFoundError:
    formula_df = get_formula_df(token, PIMP_HOST, analysis_id, database_name, polarity='positive')
    formula_df.to_pickle(formula_df_filename)

experimental_design_filename = os.path.join(os.getcwd(), 'test_data', 'experimental_design.p')
try:
    experimental_design_filename = os.path.join(os.getcwd(), 'test_data', 'experimental_design.p')
    with open(experimental_design_filename, 'rb') as f:
        experimental_design = pickle.load(f)
except FileNotFoundError:
    experimental_design = get_experimental_design(token, PIMP_HOST, analysis_id)
    with open(experimental_design_filename, 'wb') as f:
        pickle.dump(experimental_design, f)

ds = DataSource(int_df, formula_df, experimental_design, DATABASE_PIMP_KEGG)
pals = PALS(ds, min_replace=5000)
pathway_df = pals.get_pathway_df()
output = os.path.join(os.getcwd(), 'test_data', 'pathway_df_pimp_kegg.csv')
pathway_df.to_csv(output)