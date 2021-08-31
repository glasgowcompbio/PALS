#!/usr/bin/env python
import getopt
import sys
import time
from argparse import ArgumentParser
from pathlib import Path

from mummichog.reporting import LocalExporting

sys.path.append('.')

from pals.ORA import ORA
from pals.PLAGE import PLAGE
from pals.GSEA import GSEA
from pals.common import *
from pals.feature_extraction import DataSource
from mummichog.get_user_data import *
from mummichog.functional_analysis import *
import pandas.io.common

def pair(arg):
    # https://stackoverflow.com/questions/33499173/how-to-create-argument-of-type-list-of-pairs-with-argparse
    # For simplity, assume arg is a pair of integers
    # separated by a slash. If you want to do more
    # validation, raise argparse.ArgumentError if you
    # encounter a problem.
    return [str(x) for x in arg.split('/')]


def cli_options(opts):
    '''
    Ongoing work in version 2, making some options obsolete.

    obsolete parameters:
    'analysis': 'total',
    'targeted': False,
    'evidence': 3,
    'visualization': 2,

    '''
    time_stamp = str(time.time())

    optdict = {
        'cutoff': 0,

        'method': '',
        'network': 'human_mfn',
        'modeling': None,

        'mode': 'pos_default',
        'instrument': 'unspecified',
        'force_primary_ion': True,

        'workdir': '',
        'input': '',
        'reference': '',
        'infile': '',
        'output': '',
        'permutation': 100,
        'outdir': 'mcgresult' + time_stamp,


        #PALS OPTIONS

        'database': DATABASE_REACTOME_KEGG,
        'comparisons': '',
        'min replace': SMALL,
        'min hits': MIN_HITS,
        'species': 'Homo sapiens',
        'use all reactome pathways': False,
        'connect to reactome server': False

    }
    booleandict = {'T': True, 'F': False, 1: True, 0: False,
                   'True': True, 'False': False, 'TRUE': True, 'FALSE': False, 'true': True, 'false': False,
                   }
    modedict = {'default': 'pos_default', 'pos': 'pos_default', 'pos_default': 'pos_default',
                'dpj': 'dpj_positive', 'positive': 'generic_positive', 'Positive': 'generic_positive',
                'negative': 'negative', 'Negative': 'negative',
                }
    # update default from user argument
    for o, a in opts:
        if o in ("-x", "--method"):
            optdict['method'] = a
        elif o in ("-a", "--analysis"):
            optdict['analysis'] = a
        elif o in ("-c", "--cutoff"):
            optdict['cutoff'] = float(a)
        elif o in ("-t", "--targeted"):
            optdict['targeted'] = booleandict.get(a, False)
        elif o in ("-n", "--network"):
            optdict['network'] = a
        elif o in ("-z", "--force_primary_ion"):
            optdict['force_primary_ion'] = booleandict.get(a, True)
        elif o in ("-d", "--modeling"):
            optdict['modeling'] = a
        elif o in ("-e", "--evidence"):
            optdict['evidence'] = int(a)
        elif o in ("-m", "--mode"):
            optdict['mode'] = modedict.get(a, a)
        elif o in ("-u", "--instrument"):
            optdict['instrument'] = a
        elif o in ("-v", "--visualization"):
            optdict['visualization'] = int(a)
        elif o in ("-k", "--workdir"):
            optdict['workdir'] = a
        elif o in ("-i", "--input"):
            optdict['input'] = a
        elif o in ("-r", "--reference"):
            optdict['reference'] = a
        elif o in ("-f", "--infile"):
            optdict['infile'] = a
        elif o in ("-o", "--output"):
            optdict['output'] = a.replace('.csv', '')
            optdict['outdir'] = '.'.join([time_stamp, a.replace('.csv', '')])

        elif o in ("-p", "--permutation"):
            optdict['permutation'] = int(a)


        #PALS OPTIONS

        elif o in ("-q", "--comparisons"):
            optdict['comparisons'] = a
        elif o in ("-s", "--database_name"):
            optdict['database'] = a
        elif o in ("-j", "--min_replace"):
            optdict['min replace'] = a
        elif o in ('-b', '--min_hits'):
            optdict['min hits'] = a
        elif o in ("-w", "--species"):
            optdict['species'] = a
        elif o in ('-y', '--use_all_reactome_pathways'):
            optdict['use all reactome pathways'] = a
        elif o in ('-g', '--connect_to_reactome'):
            optdict['connect to reactome server'] = a

        else:
            print("Unsupported argument ", o)

    return optdict


def dispatcher():
    '''
    Dispatch command line arguments to corresponding functions.
    No user supplied id is used in version 1.
    User supplied IDs, str_mz_rtime IDs and targeted metabolites will be supported in version 2.

    '''
    helpstr = '''
    Usage example:
    
    -x --method {MUMMICHOG, PLAGE, ORA, GSEA)
    
    MUMMICHOG:
    
    python main.py -f mydata.txt -o myoutput
    
    
        -f, --infile: single file as input, 
              containing all features with tab-delimited columns
              m/z, retention time, p-value, statistic score

        -n, --network: network model to use (default human_mfn; models being ported to version 2), 
              [human_mfn, worm]

        -o, --output: output file identification string (default 'mcgresult')
        -k, --workdir: directory for all data files.
              Default is current directory.

        -m, --mode: analytical mode of mass spec, [positive, negative, pos_defult].
              Default is pos_defult, a short version of positive.
        -u, --instrument: Any integer, treated as ppm of instrument accuracy. Default is 10. 

        -p, --permutation: number of permutation to estimate null distributions.
              Default is 100.
        -z,   --force_primary_ion: one of primary ions, 
              ['M+H[1+]', 'M+Na[1+]', 'M-H2O+H[1+]', 'M-H[-]', 'M-2H[2-]', 'M-H2O-H[-]'],  
              must be present for a predicted metabolite, [True, False].
              Default is True.

        -c, --cutoff: optional cutoff p-value in user supplied statistics,
              used to select significant list of features. 
        -d, --modeling: modeling permutation data, [no, gamma].
              Default is no.
              
    PALS (ORA, PLAGE, MSEA):
    
        POSITONAL ARGUMENTS:
        [-h] intensity_csv annotation_csv output file comparisons
        
        FLAGS:
        -s --database_name: name of source database 
             {PiMP_KEGG,COMPOUND,ChEBI,UniProt,ENSEMBL}
        
        -j --min_replace: minimum intensity value for data imputation. 
             Defaults to 5000
        
        -w --species: species name for reactome pathway query e.g --species "Homo sapiens". 
             Defaults to Homo sapiens
        
        -y --use_all_reactome_pathways: option to use all pathways for Reactome pathway query. 
             If this is not used only metabolic pathways will be queried
        
        -g --connect_to_reactome_server: Whether to connect to an instance of Neo4j server hosting Reactome database. 
             If not specified then offline mode (using a downloaded copy of selected Reactome pathways will be used
    
        '''

    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:c:t:d:e:m:n:u:z:v:k:i:r:f:o:p:x:q:s:j:b:w:y:g:",
                                   ["analysis=", "cutoff", "targeted=", "modeling=", "evidence=", "mode=",
                                    "network=", "instrument=", "force_primary_ion",
                                    "visualization=", "workdir=", "input=",
                                    "reference=", "infile=", "output=", "permutation=", 'method', 'comparisons=',
                                    'database', 'min replace', 'min hits', 'species', 'use all reactome pathways', 'connect to reactome server'])
        if not opts:
            print(helpstr)
            sys.exit(2)

    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    return cli_options(opts)


def main():

    option_dictionary = dispatcher()
    method = None

    if option_dictionary['method'] == 'mummichog':
        user_data = InputUserData(option_dictionary)

        if user_data.paradict['network'] in ['human', 'hsa', 'Human', 'human_mfn', 'hsa_mfn', '']:
            theoretical_model = metabolicNetwork(metabolicModels[ 'human_model_mfn' ])

        elif user_data.paradict['network'] in ['worm', 'C. elegans', 'icel1273', 'Caenorhabditis elegans']:
            theoretical_model = metabolicNetwork(metabolicModels[ 'worm_model_icel1273' ])

        else:
            raise KeyError ("Unsupported species and/or model")

        mixed_network = DataMeetModel(theoretical_model, user_data)

        PA = PathwayAnalysis(mixed_network.model.metabolic_pathways, mixed_network)
        PA.cpd_enrich_test()

        MA = ModularAnalysis(mixed_network)
        MA.dispatch()

        AN = ActivityNetwork(mixed_network, set(PA.collect_hit_Trios() + MA.collect_hit_Trios()))

        Local = LocalExporting(mixed_network, PA, MA, AN)
        Local.run()



    elif option_dictionary['method'].lower() == 'plage' or option_dictionary['method'].lower() == 'ora' or option_dictionary['method'].lower() == 'gsea':

        intensity_csv = sys.argv[3]
        annotation_csv = sys.argv[4]
        output = sys.argv[5]
        comparisons = sys.argv[6:8]

        print()

        int_df, annotation_df, groups = load_data(intensity_csv, annotation_csv)

        experimental_design = {
            'groups': groups,
            'comparisons': []
        }

        for case, control in sys.argv[6:8]:
            experimental_design['comparisons'].append({
                'case': case,
                'control': control,
                'name': '%s/%s' % (case, control)
            })

        # fill in other variables

        database_name = option_dictionary['database']
        min_replace = option_dictionary['min replace']
        min_hits = option_dictionary['min hits']
        reactome_species = option_dictionary['species']
        reactome_metabolic_pathway_only = not option_dictionary['use all reactome pathways']
        reactome_query = option_dictionary['connect to reactome server']

        ds = DataSource(int_df, annotation_df, experimental_design, database_name,
                        reactome_species=reactome_species,
                        reactome_metabolic_pathway_only=reactome_metabolic_pathway_only,
                        reactome_query=reactome_query, min_replace=min_replace, min_hits=min_hits)

        # run the selected pathway analysis method
        method = None
        if option_dictionary['method'].upper() == PATHWAY_ANALYSIS_PLAGE:
            method = PLAGE(ds)
        elif option_dictionary['method'].upper() == PATHWAY_ANALYSIS_ORA:
            method = ORA(ds)
        elif option_dictionary['method'].upper() == PATHWAY_ANALYSIS_GSEA:
            method = GSEA(ds)
        assert method is not None

    # save the results
    #df = method.get_results()

    # filter results to show only the columns we want
    #try:
        #df.drop(columns=['sf', 'exp_F', 'Ex_Cov'], inplace=True)
    #except KeyError:
        #pass

    # drop all columns containing 'comb_p' since now they should be the same as the uncombined columns
    # (after we have removed the hypergeometric weights)
    #df = df[df.columns.drop(list(df.filter(regex='comb_p')))]

    # sort by pw_name case-insensitive
    # https://stackoverflow.com/questions/41656623/pandas-dataframe-sort-ignoring-the-case
    #df['temp'] = df['pw_name'].str.upper()
    #df.sort_values('temp', inplace=True)
    #del df['temp']


    #logger.info('Saving PALS results to %s' % output)
    #df.to_csv(output)


if __name__ == '__main__':
    main()
