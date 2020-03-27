#!/usr/bin/env python

import sys
from urllib.parse import quote

import numpy as np
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
    st.write('To begin, please upload your intensity ([example](https://raw.githubusercontent.com/glasgowcompbio/PALS/'
             'master/notebooks/test_data/beer/int_df.csv)) and annotation ([example](https://raw.githubusercontent.com/'
             'glasgowcompbio/PALS/master/notebooks/test_data/beer/annotation_df.csv)) matrices from the sidebar. Next, '
             'select a case and control group, the pathway analysis method as well as the database to use. PALS will '
             'perform pathway analysis and display the results below.')

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
        ds = get_data_source(annotation_df, database_name, experimental_design, int_df, min_replace,
                             reactome_metabolic_pathway_only, reactome_query, reactome_species)

        if len(ds.dataset_pathways) == 0:
            st.error('No matching pathways found for this data. Ensure that the database and species are correct.')
        else:
            df = pathway_analysis(significant_column, selected_method, ds)
            st.success('Pathway analysis successful!')

            token = None
            if use_reactome:
                token = send_expression_data(ds, case, control, reactome_species)
            process_results(significant_column, df, use_reactome, token)


@st.cache
def get_data_source(annotation_df, database_name, experimental_design, int_df, min_replace,
                    reactome_metabolic_pathway_only, reactome_query, reactome_species):
    ds = DataSource(int_df, annotation_df, experimental_design, database_name,
                    reactome_species=reactome_species,
                    reactome_metabolic_pathway_only=reactome_metabolic_pathway_only,
                    reactome_query=reactome_query, min_replace=min_replace)
    return ds


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


def process_results(significant_column, df, use_reactome, token):
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
                               step=0.05)
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
        st.header('Pathway Browser')
        st.write('To display additional information on significantly changing pathways, please select them in'
                 ' the list below. Entries in the list are in ascending order of their pathway activity p-values.')

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

            status_code, json_response = get_reactome_info(stId)
            if status_code == 200:
                # logger.debug(json_response)

                # st.subheader(pw)
                label = '%s: %s' % (stId, pw_name)
                info_url = 'https://reactome.org/content/detail/%s' % stId
                header_markdown = '### %s [[info]](%s)' % (label, info_url)
                if token is not None:
                    viewer_url = 'https://reactome.org/PathwayBrowser/#/%s&DTAB=AN&ANALYSIS=%s' % (stId, token)
                    header_markdown += ' [[viewer]](%s)' % viewer_url
                st.write(header_markdown)

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

                image_url = 'https://reactome.org/ContentService/exporter/diagram/%s.png?quality=8' \
                            '&diagramProfile=standard&analysisProfile=strosobar' % stId
                if token is not None:
                    image_url += '&token=%s&resource=TOTAL&expColumn=0' % token
                logger.debug('image_url = %s' % image_url)
                st.image(image_url, use_column_width=True)

                # print reactions
                # for event in json_response['hasEvent']:
                #     name = event['name'][0]
                #     species = event['speciesName']
                #     event_str = '- %s (%s)' % (name, species)
                #     st.write(event_str)


def send_expression_data(ds, case, control, species):
    # send expression data to reactome for diagram exporter
    int_df = ds.change_zero_peak_ints()
    annot_df = ds.get_annotations()
    design = ds.get_experimental_design()

    case_cols = design['groups'][case]
    control_cols = design['groups'][control]

    df = annot_df.join(int_df).set_index('entity_id').sort_index()
    case_df = np.log2(df[case_cols])
    control_df = np.log2(df[control_cols])
    lfcs = np.mean(case_df, axis=1) - np.mean(control_df, axis=1)

    # if there are multiple values for each compound, take the largest fold-change (in either direction) to display
    # TODO: we need a better way to do this
    fold_changes = {}
    for idx, lfc in lfcs.iteritems():
        if idx not in fold_changes:
            fold_changes[idx] = lfc
        else:
            if abs(lfc) > abs(fold_changes[idx]):
                fold_changes[idx] = lfc

    expression_df = pd.DataFrame.from_dict(fold_changes.items()).set_index(0)
    expression_df = expression_df.rename(columns={1: 'log_FC'})
    expression_df.index.name = '#id'

    expression_data = expression_df.to_csv(sep='\t', header=True, index_label='#id', float_format='%.15f')
    encoded_species = quote(species)

    status_code, json_response = send_reactome_expression_data(expression_data, encoded_species)
    if status_code == 200:
        pathways_df, reactome_url, reactome_token = parse_reactome_json(json_response)
        return reactome_token
    else:
        st.warning('Failed to submit expression data to Reactome.org (status_code=%d)' % status_code)
        return None


@st.cache
def get_reactome_info(stId):
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


@st.cache
def send_reactome_expression_data(data, encoded_species):
    # refer to https://reactome.org/AnalysisService/#/identifiers/getPostTextUsingPOST
    url = 'https://reactome.org/AnalysisService/identifiers/?interactors=false&species=' + encoded_species + \
          '&sortBy=ENTITIES_PVALUE&order=ASC&resource=TOTAL&pValue=1&includeDisease=true'
    logger.debug('POSTing expression data to Reactome Analysis Service: ' + url)

    # make a POSt request to Reactome Analysis service
    response = requests.post(url, headers={'Content-Type': 'text/plain'}, data=data.encode('utf-8'))
    logger.debug('Received HTTP status code: %d' % response.status_code)

    status_code = response.status_code
    if status_code == 200:
        json_response = json.loads(response.text)
    else:
        json_response = None
    return status_code, json_response


def parse_reactome_json(json_response):
    # see https://reactome.org/userguide/analysis for results explanation
    token = json_response['summary']['token']
    pathways = json_response['pathways']

    reactome_url = 'https://reactome.org/PathwayBrowser/#DTAB=AN&ANALYSIS=' + token
    logger.debug('Received expression analysis token: ' + token)

    # https://stackoverflow.com/questions/6027558/flatten-nested-dictionaries-compressing-keys
    pathways_df = pd.json_normalize(pathways, sep='_')
    return pathways_df, reactome_url, token


max_width()
main()