#!/usr/bin/env python
import sys

import streamlit as st

sys.path.append('.')

from pals.run_gui_pathway import show_pathway_widgets, run_pathway_analysis, process_results, show_results
from pals.common import *

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
    # st.write(
    #     'Understanding changing pathways or molecular families can be incredibly useful in the interpretation and '
    #     'understanding of complex datasets from metabolomics experiments. PALS is a Python package to perform the '
    #     'ranking of significantly-changing metabolite set in different experimental conditions through the decomposition '
    #     'of activity levels calculated from peak intensities.')
    st.write(
        'Please select an analysis type from the sidebar. Two types of analysis are available, depending on how '
        'metabolites are grouped into sets. You can perform either an analysis of **pathways** (metabolite sets '
        'grouped by chemical reactions) or of **molecular families** (metabolite sets grouped their fragmentation '
        'spectral). For pathway analysis, PALS Viewer supports KEGG and Reactome databases. For molecular family '
        'analysis, PALS Viewer supports the output from Feature-based Molecular Network (FBMN) from GNPS.'
    )

    image_url = 'https://raw.githubusercontent.com/glasgowcompbio/PALS/master/images/logo_transparent.png'
    st.sidebar.image(image_url, use_column_width=True)
    st.sidebar.header('Analysis Type')
    analysis_type = st.sidebar.radio(
        "Please choose an analysis to perform",
        (GUI_PATHWAY_ANALYSIS, GUI_MOLECULAR_FAMILY_ANALYSIS))

    if analysis_type == GUI_PATHWAY_ANALYSIS:

        st.subheader(GUI_PATHWAY_ANALYSIS)
        st.write('Please upload your intensity ([example](https://github.com/glasgowcompbio/PALS/raw/master/'
                 'notebooks/test_data/HAT/int_df.csv)) and annotation ([example](https://raw.githubusercontent.com/'
                 'glasgowcompbio/PALS/master/notebooks/test_data/HAT/annotation_df.csv)) matrices from the sidebar. Next, '
                 'select a case and control group, the pathway analysis method as well as the database to use. '
                 'Click the **Run** button. PALS will perform pathway analysis and display the results below.')

        st.sidebar.subheader('Input Files')
        intensity_csv = st.sidebar.file_uploader("Choose an intensity CSV file", type=['txt', 'csv'])
        annotation_csv = st.sidebar.file_uploader("Choose an annotation CSV file", type=['txt', 'csv'])
        if intensity_csv is not None and annotation_csv is not None:  # data is loaded
            params = show_pathway_widgets(intensity_csv, annotation_csv)
            significant_column = '%s/%s p-value' % (params['case'], params['control'])

            # run a pathway analysis method, e.g. PLAGE
            df, token = run_pathway_analysis(params)

            # display pathway analysis results when available
            if df is not None:
                df = process_results(df, significant_column)
                show_results(df, params['use_reactome'], token)


max_width()
main()
