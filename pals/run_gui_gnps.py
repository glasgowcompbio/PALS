import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st

from pals.PLAGE import PLAGE
from pals.common import *
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


@st.cache
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

    all_groups, all_samples, entity_dict, intensities_df, dataset_pathways_to_row_ids = get_plot_data(ds)
    results = {
        'df': df,
        'all_groups': all_groups,
        'all_samples': all_samples,
        'entity_dict': entity_dict,
        'intensities_df': intensities_df,
        'dataset_pathways_to_row_ids': dataset_pathways_to_row_ids
    }
    return results


@st.cache(suppress_st_warning=True)
def fetch_GNPS_data(gnps_url, metadata_df, comparisons):
    database_name = DATABASE_GNPS
    loader = GNPSLoader(database_name, gnps_url, metadata_df, comparisons)
    database = loader.load_data()
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
    my_bar.empty()
    return df


@st.cache
def get_plot_data(gnps_ds):
    experimental_design = gnps_ds.get_experimental_design()
    all_samples = []
    all_groups = []
    for group in experimental_design['groups']:
        samples = experimental_design['groups'][group]
        all_samples.extend(samples)
        all_groups.extend([group] * len(samples))
    entity_dict = gnps_ds.entity_dict
    intensities_df = gnps_ds.standardize_intensity_df()
    dataset_pathways_to_row_ids = gnps_ds.dataset_pathways_to_row_ids
    return all_groups, all_samples, entity_dict, intensities_df, dataset_pathways_to_row_ids


@st.cache
def process_gnps_results(df, significant_column):
    # filter results to show only the columns we want
    try:
        df = df.drop(columns=['sf', 'exp_F', 'Ex_Cov', 'unq_pw_F', 'F_coverage'])
    except KeyError:
        pass
    df = df[df.columns.drop(list(df.filter(regex='comb_p')))]

    # sort column
    count_col = 'tot_ds_F'
    df = df.sort_values([significant_column, count_col], ascending=[True, False])

    # reorder and rename columns
    df = df[['pw_name', significant_column, 'tot_ds_F']]

    df = df.rename(columns={
        'pw_name': 'Components',
        significant_column: 'p-value',
        'tot_ds_F': 'No. of members',
    })
    return df


def show_gnps_results(df, results):
    # unpack results
    all_groups = results['all_groups']
    all_samples = results['all_samples']
    entity_dict = results['entity_dict']
    intensities_df = results['intensities_df']
    dataset_pathways_to_row_ids = results['dataset_pathways_to_row_ids']

    # write header -- metabolite family ranking
    st.header('Molecular Family Ranking')
    st.write(
        ' The following table shows a ranking of molecular families ("components") based on their activity levels. '
        'Entries in the table can be filtered by p-values and the number of members ("clusters") for each '
        'molecular family.')

    # filter by significant p-values
    pval_threshold = st.slider('Filter molecular families with p-values less than', min_value=0.0, max_value=1.0,
                               value=0.05,
                               step=0.05)
    df = df[df['p-value'] <= pval_threshold].copy()

    # filter by formula hits
    min_hits = 1
    max_hits = max(df['No. of members'])
    formula_threshold = st.slider('Filter molecular families having members at least', min_value=min_hits,
                                  max_value=max_hits,
                                  value=10, step=1)
    df = df[df['No. of members'] >= formula_threshold]

    st.markdown(get_table_download_link(df), unsafe_allow_html=True)
    st.write(df)

    # write header -- pathway info
    st.header('Molecular Family Browser')
    st.write('To display additional information on significantly changing molecular families, please select them in'
             ' the list below. Entries are listed in ascending order according to their activity p-values and the'
             ' number of members.')

    choices = []
    for idx, row in df.iterrows():
        pw_name = row['Components']
        p_value = row['p-value']
        no_members = row['No. of members']
        choices.append('%s (p-value=%.4f, members=%d)' % (pw_name, p_value, no_members))
    selected = st.selectbox(
        'Select molecular family', choices)

    tokens = selected.split(' ')
    idx = tokens[2][1:]
    row = df.loc[idx, :]
    # st.write(row)

    members = dataset_pathways_to_row_ids[idx]
    member_df = get_member_df(entity_dict, members)
    plot_heatmap(all_groups, all_samples, intensities_df, member_df, members, row)
    display_member_df(member_df)


def get_member_df(entity_dict, members):
    # get group info
    # print('%s p-value=%.4f' % (pw_name, p_value))
    data = []
    for member in members:
        member_info = entity_dict[member]
        unique_id = member_info['unique_id']
        library_id = member_info['LibraryID']
        gnps_linkout_network = member_info['GNPSLinkout_Network']
        no_spectra = member_info['number of spectra']
        rt = member_info['RTConsensus']
        mz = member_info['precursor mass']
        intensity = member_info['SumPeakIntensity']
        temp = [unique_id, library_id, mz, rt, intensity, no_spectra, gnps_linkout_network]
        data.append(temp)
    member_df = pd.DataFrame(data, columns=['id', 'LibraryID', 'Precursor m/z', 'RTConsensus', 'PrecursorInt',
                                            'no_spectra', 'link']).set_index('id')
    return member_df


def display_member_df(member_df):
    st.subheader('Members')
    member_df['link'] = member_df['link'].apply(make_clickable)
    # https://discuss.streamlit.io/t/display-urls-in-dataframe-column-as-a-clickable-hyperlink/743/5
    st.write(member_df.to_html(escape=False), unsafe_allow_html=True)


def plot_heatmap(all_groups, all_samples, intensities_df, member_df, members, row):
    # Create a categorical palette to identify the networks
    used_groups = list(set(all_groups))
    group_pal = sns.husl_palette(len(used_groups), s=.45)
    group_lut = dict(zip(map(str, used_groups), group_pal))

    # Convert the palette to vectors that will be drawn on the side of the matrix
    group_intensities = intensities_df.loc[members][all_samples]
    group_colours = pd.Series(all_groups, index=group_intensities.columns).map(group_lut)
    group_colours.name = 'groups'

    # plot heatmap
    g = sns.clustermap(group_intensities, center=0, cmap='vlag', col_colors=group_colours,
                       col_cluster=False, linewidths=0.75, cbar_pos=(1.0, 0.3, 0.05, 0.5))
    pw_name = row['Components']
    plt.suptitle('%s' % (pw_name), fontsize=24, y=0.9)

    # draw group legend
    for group in used_groups:
        g.ax_col_dendrogram.bar(0, 0, color=group_lut[group], label=group, linewidth=0)
    g.ax_col_dendrogram.legend(loc="right")

    # make the annotated peaks to have labels in bold
    annotated_df = member_df[member_df['LibraryID'].notnull()]
    annotated_peaks = annotated_df.index.values
    for label in g.ax_heatmap.get_yticklabels():
        if label.get_text() in annotated_peaks:
            label.set_weight("bold")
            label.set_color("green")
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # render plot
    st.pyplot()


def make_clickable(link):
    # target _blank to open new window
    # extract clickable text to display for your link
    text = 'GNPSLinkout_Network'
    return f'<a target="_blank" href="{link}">{text}</a>'
