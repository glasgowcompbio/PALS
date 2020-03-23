#!/usr/bin/env python

import sys

import streamlit as st

sys.path.append('.')

from pals.ORA import ORA
from pals.PALS import PALS
from pals.common import *
from pals.feature_extraction import DataSource


def main():
    st.title('Pathway Activity Level Scoring (PALS)')

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
            (DATABASE_PIMP_KEGG, DATABASE_REACTOME_KEGG, DATABASE_REACTOME_CHEBI, DATABASE_REACTOME_UNIPROT,
             DATABASE_REACTOME_ENSEMBL)
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
            reactome_query = st.sidebar.checkbox('Connect to a Neo4j (Reactome) database.', value=False)

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
            process_results(significant_column, df)


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


def process_results(significant_column, df):
    # filter by significant p-values
    threshold = st.slider('Specify p-value threshold', min_value=0.0, max_value=1.0, value=0.05, step=0.01)
    df = df[df[significant_column] <= threshold]
    st.write(df)


main()
