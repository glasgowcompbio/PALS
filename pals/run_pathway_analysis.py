#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.append('.')

from pals.ORA import ORA
from pals.PALS import PALS
from pals.common import *
from pals.feature_extraction import DataSource


def pair(arg):
    # https://stackoverflow.com/questions/33499173/how-to-create-argument-of-type-list-of-pairs-with-argparse
    # For simplity, assume arg is a pair of integers
    # separated by a slash. If you want to do more
    # validation, raise argparse.ArgumentError if you
    # encounter a problem.
    return [str(x) for x in arg.split('/')]


def build_parser():
    # example usage:
    # run_pals.py PALS intensity.csv annotation.csv --comparisons beer1/beer2 beer3/beer4 --db DATABASE_PIMP_KEGG --min_replace 5000 --plage_weight 5 --hg_weight 1
    # run_pals.py ORA intensity.csv annotation.csv --comparisons beer1/beer2 beer3/beer4 --db DATABASE_REACTOME_KEGG --species Homo sapiens --metabolic-pathways-only True --reactome-query True

    parser = ArgumentParser(description="Run Pathway Analysis")

    # required parameters
    parser.add_argument('method', default=PATHWAY_ANALYSIS_PALS, help='Pathway Analysis Method',
                        choices=(PATHWAY_ANALYSIS_PALS, PATHWAY_ANALYSIS_ORA)),
    parser.add_argument('intensity_csv', type=Path, help='Intensity CSV file')
    parser.add_argument('annotation_csv', type=Path, help='Annotation CSV file')
    parser.add_argument('output_file', type=Path, help='PALS analysis output')
    parser.add_argument('--db', required=True, default=DATABASE_REACTOME_KEGG, help='Database name', choices=(
        DATABASE_PIMP_KEGG, DATABASE_REACTOME_KEGG, DATABASE_REACTOME_CHEBI, DATABASE_REACTOME_UNIPROT,
        DATABASE_REACTOME_ENSEMBL))
    parser.add_argument('--comparisons', required=True, type=pair, nargs='+')

    # common parameters
    parser.add_argument('--min_replace', type=float, default=5000.0,
                        help='Minimum intensity of MS1 peaks for data imputation  (default: %(default)s)')
    parser.add_argument('--plage_weight', type=float, default=5.0, help=' (default: %(default)s)')
    parser.add_argument('--hg_weight', type=float, default=1.0, help=' (default: %(default)s)')

    # reactome parameters
    parser.add_argument('--species', default=REACTOME_SPECIES_HOMO_SAPIENS, help='Species name',
                        choices=(
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
                        ))
    parser.add_argument('--use_all_reactome_pathways', default=False, action='store_true',
                        help='Limit pathway queries to Reactome metabolic pathways only.')
    parser.add_argument('--connect_to_reactome_server', default=False, action='store_true',
                        help='Query pathways by connecting to a Neo4j server hosting Reactome database.')
    return parser


def load_data(intensity_csv, annotation_csv):
    """
    Loads PALS data from csv files
    :param intensity_csv: a CSV file of peak intensities
    :param annotation_csv: a CSV file of peak annotation
    :return: the peak intensity df, annotation df, and also grouping information
    """
    # load intensity dataframe from csv
    # skip the second row containing group information
    # set the first column (peak ids) in csv to be the index
    int_df = pd.read_csv(intensity_csv, skiprows=[1], index_col=0)
    logger.debug('Loaded %d x %d peak intensities from %s' % (int_df.shape[0], int_df.shape[1], intensity_csv))

    # load annotation dataframe from csv, setting the first column to be the index
    annotation_df = pd.read_csv(annotation_csv, index_col=0)
    logger.debug('Loaded %d peak annotations from %s' % (annotation_df.shape[0], annotation_csv))

    # extract grouping information from csv
    groups = get_groups(intensity_csv)
    logger.debug('Loaded groups: %s' % groups)

    return int_df, annotation_df, groups


def get_groups(csv_file):
    """
    Extract group information from peak intensity csv file
    :param csv_file: An input CSV file.
    The first row in the CSV file should contain samples information.
    The second row in the CSV file should contain grouping information.
    :return: a dictionary where key is the group name and values are the samples under that group
    """
    grouping = defaultdict(list)
    with open(csv_file, 'r') as f:
        # extract the first and second lines containing samples and grouping information
        samples = None
        groups = None
        for i, line in enumerate(f):
            if i == 0:
                samples = line.strip().split(',')
            if i == 1:
                groups = line.strip().split(',')
                break

        assert samples is not None, 'Missing samples line'
        assert groups is not None, 'Missing groups line'
        for sample, group in zip(samples, groups):
            if group == 'group':  # skip the column header containing the word 'group'
                continue
            grouping[group].append(sample)
    return dict(grouping)


def main(args):
    logger.info('Welcome to Pathway Activity Level Scoring (PALS)')
    logger.debug('PARAMETERS:')
    for k, v in args.items():
        logger.debug('\t- %s: %s' % (k, v))

    # load intensity df, annotation df and group information from CSV files
    intensity_csv = args['intensity_csv']
    annotation_csv = args['annotation_csv']
    int_df, annotation_df, groups = load_data(intensity_csv, annotation_csv)

    # populate comparisons from user-provided arguments
    experimental_design = {
        'groups': groups,
        'comparisons': []
    }
    for case, control in args['comparisons']:
        experimental_design['comparisons'].append({
            'case': case,
            'control': control,
            'name': '%s/%s' % (case, control)
        })

    # extract other args
    database_name = args['db']
    min_replace = args['min_replace']
    plage_weight = args['plage_weight']
    hg_weight = args['hg_weight']
    reactome_species = args['species']
    reactome_metabolic_pathway_only = not args['use_all_reactome_pathways']
    reactome_query = args['connect_to_reactome_server']

    # create Data Source
    ds = DataSource(int_df, annotation_df, experimental_design, database_name,
                    reactome_species=reactome_species,
                    reactome_metabolic_pathway_only=reactome_metabolic_pathway_only,
                    reactome_query=reactome_query, min_replace=min_replace)

    # run the selected pathway analysis method
    method = None
    if args['method'] == PATHWAY_ANALYSIS_PALS:
        method = PALS(ds, plage_weight=plage_weight, hg_weight=hg_weight)
    elif args['method'] == PATHWAY_ANALYSIS_ORA:
        method = ORA(ds)
    assert method is not None

    # save the results
    pathway_df = method.get_pathway_df()
    output_file = args['output_file']
    logger.info('Saving PALS results to %s' % output_file)
    pathway_df.to_csv(output_file)


if __name__ == '__main__':
    parser = build_parser()
    args = parser.parse_args()
    fargs = vars(args)
    main(fargs)
