#!/usr/bin/env python

import sys

import requests
import streamlit as st

sys.path.append('.')

from pals.ORA import ORA
from pals.PALS import PALS
from pals.common import *
from pals.feature_extraction import DataSource


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
    st.title('Pathway Activity Level Scoring (PALS) :twisted_rightwards_arrows:')
    st.write(
        'Understanding changing pathways can be incredibly useful in the interpretation and understanding of complex '
        'datasets from metabolomics experiments. PALS is a Python package to perform the ranking of '
        'significantly-changing metabolite pathways in different experimental conditions through the decomposition of '
        'pathway activity levels calculated from peak intensities.')
    st.write('To begin, please upload your intensity and annotation matrices from the sidebar. Select a case and'
             'control group, the pathway analysis method as well as the database to use.')

    intensity_csv = st.sidebar.file_uploader("Choose an intensity CSV file", type=['txt', 'csv'])
    annotation_csv = st.sidebar.file_uploader("Choose an annotation CSV file", type=['txt', 'csv'])
    if intensity_csv is not None and annotation_csv is not None:  # data is loaded
        int_df, annotation_df, groups = load_data(intensity_csv, annotation_csv, gui=True)
        case = st.sidebar.selectbox(
            'Case',
            list(groups.keys()),
            index=0
        )
        control = st.sidebar.selectbox(
            'Control',
            list(groups.keys()),
            index=1
        )
        experimental_design = {
            'groups': groups,
            'comparisons': [
                {
                    'case': case,
                    'control': control,
                    'name': '%s/%s' % (case, control)
                }
            ]
        }

        selected_method = st.sidebar.selectbox(
            'Pathway Analysis Method',
            (PATHWAY_ANALYSIS_PALS, PATHWAY_ANALYSIS_ORA, PATHWAY_ANALYSIS_GSEA)
        )
        database_name = st.sidebar.selectbox(
            ('Database'),
            (DATABASE_REACTOME_KEGG, DATABASE_REACTOME_CHEBI, DATABASE_PIMP_KEGG, DATABASE_REACTOME_UNIPROT,
             DATABASE_REACTOME_ENSEMBL),
        )
        min_replace = st.sidebar.number_input('Minimum intensity for data imputation', value=MIN_REPLACE)

        use_reactome = False
        if database_name in (DATABASE_REACTOME_KEGG, DATABASE_REACTOME_CHEBI, DATABASE_REACTOME_UNIPROT,
                             DATABASE_REACTOME_ENSEMBL):
            use_reactome = True

        reactome_species = None
        reactome_metabolic_pathway_only = False
        reactome_query = False
        if use_reactome:
            species_list = (
                REACTOME_SPECIES_ARABIDOPSIS_THALIANA,
                REACTOME_SPECIES_BOS_TAURUS,
                REACTOME_SPECIES_CAENORHABDITIS_ELEGANS,
                REACTOME_SPECIES_CANIS_LUPUS_FAMILIARIS,
                REACTOME_SPECIES_DANIO_RERIO,
                REACTOME_SPECIES_DICTYOSTELIUM_DISCOIDEUM,
                REACTOME_SPECIES_DROSOPHILA_MELANOGASTER,
                REACTOME_SPECIES_GALLUS_GALLUS,
                REACTOME_SPECIES_HOMO_SAPIENS,
                REACTOME_SPECIES_MUS_MUSCULUS,
                REACTOME_SPECIES_ORYZA_SATIVA,
                REACTOME_SPECIES_RATTUS_NORVEGICUS,
                REACTOME_SPECIES_SACCHAROMYCES_CEREVISIAE,
                REACTOME_SPECIES_SUS_SCROFA
            )
            human_idx = species_list.index(REACTOME_SPECIES_HOMO_SAPIENS)
            reactome_species = st.sidebar.selectbox(
                'Species',
                species_list,
                index=human_idx
            )
            reactome_metabolic_pathway_only = st.sidebar.checkbox('Limit to metabolic pathways only.', value=True)
            reactome_query = st.sidebar.checkbox('Connect to a Reactome Neo4j database (Online Mode).', value=False)

        significant_column = '%s/%s p-value' % (case, control)
        ds = DataSource(int_df, annotation_df, experimental_design, database_name,
                        reactome_species=reactome_species,
                        reactome_metabolic_pathway_only=reactome_metabolic_pathway_only,
                        reactome_query=reactome_query, min_replace=min_replace)

        if len(ds.dataset_pathways) == 0:
            st.error('No matching pathways found for this data. Ensure that the database and species are correct.')
        else:
            df = pathway_analysis(significant_column, selected_method, ds)
            st.success('Pathway analysis successful!')
            process_results(significant_column, df, use_reactome)


@st.cache
def pathway_analysis(significant_column, selected_method, ds):
    method = None
    if selected_method == PATHWAY_ANALYSIS_PALS:
        method = PALS(ds)
    elif selected_method == PATHWAY_ANALYSIS_ORA:
        method = ORA(ds)
    assert method is not None

    df = method.get_pathway_df()

    # filter results to show only the columns we want
    try:
        df = df.drop(columns=['sf', 'exp_F', 'Ex_Cov'])
    except KeyError:
        pass
    df = df[df.columns.drop(list(df.filter(regex='comb_p')))]

    # sort column
    df = df.sort_values(significant_column)
    return df


@st.cache
def send_to_reactome(stId):
    # refer to https://reactome.org/dev/content-service
    url = 'https://reactome.org/ContentService/data/query/%s' % stId
    logger.debug('Reactome URL: ' + url)

    # make a GET request to Reactome Content service
    response = requests.get(url)
    logger.debug('Response status code = %d' % response.status_code)

    status_code = response.status_code
    if status_code == 200:
        json_response = json.loads(response.text)
    else:
        json_response = None
    return status_code, json_response


def process_results(significant_column, df, use_reactome):
    # reorder and rename columns
    df = df[['pw_name', significant_column, 'tot_ds_F', 'unq_pw_F', 'F_coverage']]
    df = df.rename(columns={
        'pw_name': 'Pathways',
        'unq_pw_F': 'Pathway Formula',
        'tot_ds_F': 'Formula Hits',
        'F_coverage': 'Formula Coverage (%)'
    })

    # write header -- pathway ranking
    st.header('Pathway Ranking')
    st.write(' The following table shows a ranking of pathways based on their activity levels. Entries in the table can'
             ' be filtered by p-values and the number of formula hits for each pathway.')

    # filter by significant p-values
    pval_threshold = st.slider('Filter pathways with p-values less than', min_value=0.0, max_value=1.0, value=0.05,
                               step=0.01)
    df = df[df[significant_column] <= pval_threshold]

    # filter by formula hits
    min_hits = 1
    max_hits = max(df['Formula Hits'])
    formula_threshold = st.slider('Filter pathways with formula hits at least', min_value=min_hits, max_value=max_hits,
                                  value=2, step=1)
    df = df[df['Formula Hits'] >= formula_threshold]

    st.write(df)

    if use_reactome:
        # write header -- pathway info
        st.header('Pathway Information')
        st.write('The following shows additional information on the selected pathways.')

        choices = []
        for idx, row in df.iterrows():
            pw_name = row['Pathways']
            choices.append('%s (%s)' % (pw_name, idx))
        options = st.multiselect(
            'Select pathways', choices)

        for pw in options:
            tokens = pw.split('(')
            pw_name = tokens[0].strip()
            stId = tokens[1].strip()
            stId = stId[:-1]  # remove last ')' character from stId

            status_code, json_response = send_to_reactome(stId)
            if status_code == 200:
                logger.debug(json_response)

                # st.subheader(pw)
                link_label = '%s (%s)' % (pw_name, stId)
                link_url = 'https://reactome.org/content/detail/%s' % stId
                st.write('### [%s](%s)' % (link_label, link_url))

                row = df.loc[stId]
                pvalue = row[significant_column]
                num_hits = row['Formula Hits']
                subsubheader = '#### p-value: %.4f' % pvalue
                st.write(subsubheader)
                subsubheader = '#### Formula Hits: %d' % (num_hits)
                st.write(subsubheader)

                st.write('#### Summary:')
                for summation in json_response['summation']:
                    summary = summation['text']
                    st.write(summary)

                image_url = 'https://reactome.org/ContentService/exporter/diagram/%s.png?quality=8&' \
                            'diagramProfile=standard&analysisProfile=strosobar' % stId
                st.image(image_url, use_column_width=True)

                # for event in json_response['hasEvent']:
                #     name = event['name'][0]
                #     species = event['speciesName']
                #     event_str = '- %s (%s)' % (name, species)
                #     st.write(event_str)

max_width()
main()
