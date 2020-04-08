#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from pathlib import Path

sys.path.append('.')

from pals.ORA import ORA
from pals.PLAGE import PLAGE
from pals.GSEA import GSEA
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
    parser.add_argument('method', default=PATHWAY_ANALYSIS_PLAGE, help='Pathway Analysis Method',
                        choices=(PATHWAY_ANALYSIS_PLAGE, PATHWAY_ANALYSIS_ORA, PATHWAY_ANALYSIS_GSEA)),
    parser.add_argument('intensity_csv', type=Path, help='Intensity CSV file')
    parser.add_argument('annotation_csv', type=Path, help='Annotation CSV file')
    parser.add_argument('output_file', type=Path, help='PALS analysis output')
    parser.add_argument('--db', required=True, default=DATABASE_REACTOME_KEGG, help='Database name', choices=(
        DATABASE_PIMP_KEGG, DATABASE_REACTOME_KEGG, DATABASE_REACTOME_CHEBI, DATABASE_REACTOME_UNIPROT,
        DATABASE_REACTOME_ENSEMBL))
    parser.add_argument('--comparisons', required=True, type=pair, nargs='+')

    # common parameters
    parser.add_argument('--min_replace', type=float, default=MIN_REPLACE,
                        help='Minimum intensity of MS1 peaks for data imputation  (default: %(default)s)')

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
                        help='Limit to Reactome metabolic pathways only.')
    parser.add_argument('--connect_to_reactome_server', default=False, action='store_true',
                        help='Connect to a Neo4j (Reactome) database.')
    return parser


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
    if args['method'] == PATHWAY_ANALYSIS_PLAGE:
        method = PLAGE(ds)
    elif args['method'] == PATHWAY_ANALYSIS_ORA:
        method = ORA(ds)
    elif args['method'] == PATHWAY_ANALYSIS_GSEA:
        method = GSEA(ds)
    assert method is not None

    # save the results
    df = method.get_pathway_df()

    # filter results to show only the columns we want
    try:
        df.drop(columns=['sf', 'exp_F', 'Ex_Cov'], inplace=True)
    except KeyError:
        pass

    # drop all columns containing 'comb_p' since now they should be the same as the uncombined columns
    # (after we have removed the hypergeometric weights)
    df = df[df.columns.drop(list(df.filter(regex='comb_p')))]

    # sort by pw_name case-insensitive
    # https://stackoverflow.com/questions/41656623/pandas-dataframe-sort-ignoring-the-case
    df['temp'] = df['pw_name'].str.upper()
    df.sort_values('temp', inplace=True)
    del df['temp']

    output_file = args['output_file']
    logger.info('Saving PALS results to %s' % output_file)
    df.to_csv(output_file)


if __name__ == '__main__':
    parser = build_parser()
    args = parser.parse_args()
    fargs = vars(args)
    main(fargs)
