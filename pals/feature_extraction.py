import os
import xmltodict
import json

from loguru import logger
from collections import defaultdict


class DataSource(object):

    def __init__(self, int_df, formula_df, experimental_design, database_name='kegg'):
        # a dataframe of peak intensities, where rows = ms1_peak_id and columns = sample_name
        self.int_df = int_df

        # a dictionary specifying the experimental design
        self.experimental_design = experimental_design
        self.groups = dict(self.experimental_design['groups'].items())
        self.comparisons = self.experimental_design['comparisons']

        # load compound and pathway database information from file
        json_file = None
        if database_name == 'kegg':
            json_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data', 'kegg.json'))
            logger.debug('Loading %s' % json_file)

        with open(json_file) as f:
            data = json.load(f)

            # mapid -> pathway name
            self.pathway_dict = data['pw_dict']

            # mapid -> [ compound names ]
            self.pathway_cmpd_dict = data['pathway_cmpd_dict']

            # compound name -> formula, extracted from xml file
            self.cmpd_formula_dict = data['cmpd_formula_dict']

            # compound id -> formula, extracted from xml file
            self.cmpd_id_formula_dict = data['cmpd_id_formula_dict']

            # assert len(self.cmpd_formula_dict) == len(self.cmpd_id_formula_dict)

            # compound id -> [ mapids ], extracted from rdata file
            self.cmpd_id_pw_dict = data['cmpd_id_pw_dict']

            # map between pathway id to compound ids and formulas
            pw_cmpd_id_dict = defaultdict(list) # mapid -> [ compound ids ]
            pw_cmpd_formula_dict = defaultdict(set) # mapid -> [ formulas ]
            for cmpd_id, mapids in self.cmpd_id_pw_dict.items():
                try:
                    cmpd_formula = self.cmpd_id_formula_dict[cmpd_id]
                    for mapid in mapids:
                        pw_cmpd_id_dict[mapid].append(cmpd_id)
                        pw_cmpd_formula_dict[mapid].add(cmpd_formula)
                except KeyError:
                    continue  # TODO: skip this compound since we can't find its name or formula. Should never happen!!
            self.pw_cmpd_id_dict = dict(pw_cmpd_id_dict)
            self.pw_cmpd_formula_dict = dict(pw_cmpd_formula_dict)

            # a dataframe of peak id, originating database name, database id, formula
            formula_df = formula_df[formula_df['db'] == database_name]  # filter by db name, e.g. 'kegg'
            ds_pathways = []
            ds_pathways_peak_ids = defaultdict(list)
            for idx, row in formula_df.iterrows():
                cmpd_id = row['identifier']
                pid = idx
                try:
                    possible_pathways = self.cmpd_id_pw_dict[cmpd_id]
                    ds_pathways.extend(possible_pathways)
                    for p in possible_pathways:
                        ds_pathways_peak_ids[p].append(pid)
                except KeyError: # no information about compound -> pathway in our database
                    continue
            self.ds_pathways = set(ds_pathways)
            self.ds_pathways_peak_ids = dict(ds_pathways_peak_ids)
            self.ds_formulas = set(formula_df['formula'].values)

        # For use in the hypergeometric test - the number of unique formulas in kegg and in pathways
        # and the number of unique formulas in the ds and in pathways
        self.kegg_pw_formulas, self.unique_ds_pw_fs = self.get_unique_pw_f()

    """
    Method to return the unique pathway formulas associated with an analysis
    """

    def get_unique_pw_f(self):
        # set of compounds by id in Kegg
        kegg_cmpd_ids = set(self.cmpd_id_formula_dict.keys())
        # set of compounds by id found in the pathways
        pathway_cmpd_ids = set(self.cmpd_id_pw_dict.keys())
        # These are the pathway compound (by id) that are also present in Kegg
        pw_cmpd_ids = kegg_cmpd_ids.intersection(pathway_cmpd_ids)

        kegg_pathway_formulas = {}
        for c in pw_cmpd_ids:
            kegg_pathway_formulas[c] = self.cmpd_id_formula_dict[c]
        kegg_pw_formulas = set(kegg_pathway_formulas.values())

        # ds_formulas that are M+H or M-H - probably should use
        unique_ds_pw_fs = kegg_pw_formulas.intersection(self.ds_formulas)
        return len(kegg_pw_formulas), len(unique_ds_pw_fs)

    """ Param: the mapid of a pathway (e.g. map00010)
        Returns: the number of unique formula identifiable pathway
    """

    def get_pw_unique_F(self, mapids):
        return [len(self.pw_cmpd_formula_dict[mapid]) for mapid in mapids]

    """
    Takes in a list of pathway ids and return the number of formulas in the dataset for those pathways"""

    def get_ds_pw_compounds(self, mapids):
        num_totalF = []
        for mapid in mapids:
            if mapid in self.ds_pathways:  # get all the formulas of the pathway in the dataset
                formulas = self.pw_cmpd_formula_dict[mapid].intersection(self.ds_formulas)
                num_totalF.append(len(formulas))
            else:
                num_totalF.append(0)
        return num_totalF

    def get_comparison_samples(self, comp):
        conditions = [comp['control'], comp['case']]
        condition_samples = list(map(lambda group_name: self.groups[group_name], conditions))
        return condition_samples
