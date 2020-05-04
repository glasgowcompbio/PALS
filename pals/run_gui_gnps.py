import streamlit as st

from pals.PLAGE import PLAGE
from pals.common import *
from pals.confirm_button_hack import cache_on_button_press
from pals.feature_extraction import DataSource
from pals.loader import GNPSLoader


def show_gnps_widgets(gnps_url, metadata_csv):
    metadata_df = pd.read_csv(metadata_csv)
    choices = metadata_df['group'].unique()

    st.sidebar.subheader('Comparisons')
    control = st.sidebar.selectbox(
        'Control',
        choices,
        index=0
    )
    case = st.sidebar.selectbox(
        'Case',
        choices,
        index=1
    )
    if case == control:
        st.error("Control ('%s') cannot be the same as case ('%s')." % (control, case))
        return

    parameters = {
        'case': case,
        'control': control,
        'metadata_df': metadata_df,
        'gnps_url': gnps_url
    }
    return parameters


@cache_on_button_press('Run Analysis')
def run_gnps_analysis(params):
    case = params['case']
    control = params['control']
    comp_name = '%s/%s' % (case, control)
    comparisons = [{'case': case, 'control': control, 'name': comp_name}, ]
    gnps_url = params['gnps_url']
    metadata_df = params['metadata_df']

    # perform a POST request to get the downloadable Cytoscape results from GNPS
    # convert the retrieved data to a DataSource object
    database = fetch_GNPS_data(gnps_url, metadata_df, comparisons)

    # convert to a DataSource object that can be used by PLAGE
    ds = to_data_source(database)

    # run PLAGE decomposition on the ds
    df = PLAGE_decomposition(ds)
    p_value_col = '%s p-value' % comp_name
    count_col = 'unq_pw_F'
    df.sort_values([p_value_col, count_col], ascending=[True, False], inplace=True)
    return df


@st.cache(suppress_st_warning=True)
def fetch_GNPS_data(gnps_url, metadata_df, comparisons):
    database_name = DATABASE_GNPS
    my_bar = st.progress(0)
    loader = GNPSLoader(database_name, gnps_url, metadata_df, comparisons, streamlit_pbar=my_bar)
    database = loader.load_data()
    my_bar.progress(100)
    return database


@st.cache
def to_data_source(database):
    measurement_df = database.extra_data['measurement_df']
    annotation_df = database.extra_data['annotation_df']
    experimental_design = database.extra_data['experimental_design']
    ds = DataSource(measurement_df, annotation_df, experimental_design, None, database=database)
    return ds


@st.cache(suppress_st_warning=True)
def PLAGE_decomposition(ds):
    my_bar = st.progress(0)
    method = PLAGE(ds)
    df = method.get_pathway_df(streamlit_pbar=my_bar)
    return df
