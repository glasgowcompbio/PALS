import pandas as pd
import os

from pimp_tools import get_pimp_API_token_from_env, PIMP_HOST, get_ms1_intensities, get_ms1_peaks, get_formula_df
from feature_extraction import DataSource
from pathway_analysis import PALS

token = get_pimp_API_token_from_env()
analysis_id = 1321 # example beer analysis

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

experiment_design = {
    'groups': {
        'beer1': ['Beer_1_full1.mzXML', 'Beer_1_full2.mzXML', 'Beer_1_full3.mzXML'],
        'beer2': ['Beer_2_full1.mzXML', 'Beer_2_full2.mzXML', 'Beer_2_full3.mzXML'],
        'beer3': ['Beer_3_full1.mzXML', 'Beer_3_full2.mzXML', 'Beer_3_full3.mzXML'],
        'beer4': ['Beer_4_full1.mzXML', 'Beer_4_full2.mzXML', 'Beer_4_full3.mzXML'],
    },
    'comparisons': [
        {
            'name': 'beer1/beer2',
            'case': 'beer1',
            'control': 'beer2'
        },
        {
            'name': 'beer3/beer4',
            'case': 'beer3',
            'control': 'beer4'
        },
    ]
}

ds = DataSource(int_df, formula_df, experiment_design, database_name='kegg')
pals = PALS(ds, min_intensity=5000, num_resamples=10)
activity_df = pals.get_plage_activity_df()
plage_df = pals.set_up_resample_plage_p_df(activity_df)
pathway_df = pals.calculate_hg_values(plage_df)