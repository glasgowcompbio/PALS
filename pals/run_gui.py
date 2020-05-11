#!/usr/bin/env python
import sys

import streamlit as st

sys.path.append('.')

from pals.run_gui_pathway import show_pathway_widgets, run_pathway_analysis, process_pathway_results, \
    show_pathway_results
from pals.run_gui_gnps import show_gnps_widgets, run_analysis, process_gnps_results, show_gnps_results
from pals.common import *

# pandas display options
pd.set_option('display.max_colwidth', None)


# https://discuss.streamlit.io/t/custom-render-widths/81/8
def max_width():
    max_width_str = f"max-width: 2000px;"
    st.markdown(
        f"""
    <style>
    .reportview-container .main .block-container{{
        {max_width_str}
    }}
    </style>    
    """,
        unsafe_allow_html=True,
    )


def main():
    st.title('PALS - Pathway Activity Level Scoring')
    image_url = 'https://raw.githubusercontent.com/glasgowcompbio/PALS/master/images/logo_transparent.png'
    st.sidebar.image(image_url, use_column_width=True)
    st.sidebar.header('Analysis Type')
    analysis_type = st.sidebar.radio(
        "Please choose an analysis to perform",
        (GUI_PATHWAY_ANALYSIS, GUI_GNPS_MOLECULAR_FAMILY_ANALYSIS, GUI_GNPS_MS2LDA_ANALYSIS))

    if analysis_type == GUI_PATHWAY_ANALYSIS:

        st.sidebar.subheader('Peak Data')
        intensity_csv = st.sidebar.file_uploader("Choose an intensity CSV file", type=['txt', 'csv'])

        st.sidebar.subheader('Annotation Data')
        annotation_csv = st.sidebar.file_uploader("Choose an annotation CSV file", type=['txt', 'csv'])

        if intensity_csv is not None and annotation_csv is not None:  # data is loaded
            params = show_pathway_widgets(intensity_csv, annotation_csv)
            significant_column = '%s/%s p-value' % (params['case'], params['control'])

            # run a pathway analysis method, e.g. PLAGE
            df, token = run_pathway_analysis(params)

            # display pathway analysis results when available
            if df is not None:
                df = process_pathway_results(df, significant_column)
                show_pathway_results(df, params['use_reactome'], token)
        else:
            write_main_text(analysis_type)

    elif analysis_type == GUI_GNPS_MOLECULAR_FAMILY_ANALYSIS:

        st.sidebar.subheader('Molecular Family Data')
        gnps_url = st.sidebar.text_input('Provide a link to GNPS FBMN results')

        st.sidebar.subheader('Metadata')
        metadata_csv = st.sidebar.file_uploader("Choose a metadata CSV file", type=['txt', 'csv'])

        if len(gnps_url) > 0 and metadata_csv is not None:  # data is loaded
            params = show_gnps_widgets(gnps_url, None, metadata_csv, None)
            significant_column = '%s/%s p-value' % (params['case'], params['control'])

            # run PLAGE on the GNPS data
            results = run_analysis(params)
            df = process_gnps_results(results['df'], significant_column)
            show_gnps_results(df, results)
        else:
            write_main_text(analysis_type)

    elif analysis_type == GUI_GNPS_MS2LDA_ANALYSIS:

        st.sidebar.subheader('Motif Data')
        ms2lda_url = st.sidebar.text_input('Provide a link to GNPS-MS2LDA results')
        motif_available = len(ms2lda_url) > 0

        st.sidebar.subheader('Peak Data')
        gnps_url = st.sidebar.text_input('Provide either a link to GNPS FBMN results')
        peak_table_csv = st.sidebar.file_uploader('Or upload a peak table in CSV format', type=['txt', 'csv'])
        peak_data_available = len(gnps_url) > 0 or peak_table_csv is not None

        st.sidebar.subheader('Metadata')
        metadata_csv = st.sidebar.file_uploader("Choose a metadata CSV file", type=['txt', 'csv'])
        metadata_available = metadata_csv is not None

        if motif_available and peak_data_available and metadata_available:  # data is loaded
            params = show_gnps_widgets(gnps_url, ms2lda_url, metadata_csv, peak_table_csv)
            significant_column = '%s/%s p-value' % (params['case'], params['control'])

            # run PLAGE on the GNPS data
            results = run_analysis(params)
            df = process_gnps_results(results['df'], significant_column)
            show_gnps_results(df, results)
        else:
            write_main_text(analysis_type)


def write_main_text(analysis_type):
    st.write(
        'Please select an analysis type from the sidebar. Several types of analysis are available, depending on how '
        'metabolites are grouped into sets. For knowledge-based analysis, you can perform an analysis '
        'of **pathways** (metabolite sets grouped by chemical reactions). For analysis based on unsupervised data '
        'mining, you can perform an analysis of **molecular families** (metabolite sets clustered by fragmentation '
        'spectra) or **Mass2Motifs** (metabolite sets grouped by shared substructures). '
        'For pathway analysis, PALS supports KEGG and Reactome databases. For molecular family analysis, '
        'PALS supports the output of [Feature-based Molecular Network](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/) '
        '(FBMN) from GNPS. For Mass2Motif analysis, PALS supports the output from the '
        '[MS2LDA workflow on GNPS](https://ccms-ucsd.github.io/GNPSDocumentation/ms2lda/).'
    )

    if analysis_type == GUI_PATHWAY_ANALYSIS:
        st.subheader(GUI_PATHWAY_ANALYSIS)
        st.write('Please upload your intensity ([example](https://github.com/glasgowcompbio/PALS/raw/master/'
                 'notebooks/test_data/HAT/int_df.csv)) and annotation ([example](https://raw.githubusercontent.com/'
                 'glasgowcompbio/PALS/master/notebooks/test_data/HAT/annotation_df.csv)) CSV files from the sidebar. Next, '
                 'select the case and control groups, the pathway analysis method as well as the database to use.')

    elif analysis_type == GUI_GNPS_MOLECULAR_FAMILY_ANALYSIS:
        st.subheader(GUI_GNPS_MOLECULAR_FAMILY_ANALYSIS)
        st.write(
            'Please provide a link to your FBMN results at GNPS ([example](https://'
            'gnps.ucsd.edu/ProteoSAFe/status.jsp?task=0a8432b5891a48d7ad8459ba4a89969f)). Also provide a metadata CSV'
            ' file ([example](https://github.com/glasgowcompbio/PALS/raw/master/notebooks/test_data/AGP/'
            'AG_Plants_extremes_metadata_df.csv)) describing the assignment of samples to groups. '
            'Select the case and control groups to begin.')

    elif analysis_type == GUI_GNPS_MS2LDA_ANALYSIS:
        st.subheader(GUI_GNPS_MS2LDA_ANALYSIS)
        st.write(
            'Please provide a link to your motif-based analysis results through GNPS-MS2LDA. '
            'Next provide MS1 intensity information, either through a link to your FBMN results at GNPS or by '
            'uploading the peak table directly. Finally provide a metadata CSV file describing the assignment of '
            'samples to groups. Select the case and control groups to begin. Two example datasets are provided below.')
        st.write(
            '- American Gut Project dataset -- '
            '**motif data:** [[GNPS-MS2LDA results]](https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=7c34badae00e43bc87b195a706cf1f43), '
            '**peak data:** [[GNPS-FBMN results]](https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=0a8432b5891a48d7ad8459ba4a89969f), '
            '**metadata:** [[Metadata CSV]](https://github.com/glasgowcompbio/PALS/raw/master/notebooks/test_data/AGP/AG_Plants_extremes_metadata_df.csv)'
        )
        st.write(
            '- Rhamnaceae dataset -- '
            '**motif data:** [[GNPS-MS2LDA results]](https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=b33b2697e7924ee1920dba207ed57733), '
            '**peak data:** [[Peak table CSV]](https://github.com/glasgowcompbio/PALS/raw/master/notebooks/test_data/Rhamnaceae/171205_71extracts_MS1peaktable_MS2LDA_comma.csv), '
            '**metadata:** [[Metadata CSV]](https://github.com/glasgowcompbio/PALS/raw/master/notebooks/test_data/Rhamnaceae/MetaData_Rhamnaceae.csv)'
        )

max_width()
main()
