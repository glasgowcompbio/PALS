import os
import xmltodict
import json

from loguru import logger
from collections import defaultdict


class DataSource(object):

    def __init__(self, measurement_df, annotation_df, experimental_design, database_name='kegg'):
        # a dataframe of peak intensities, where rows = ms1_peak_id and columns = sample_name
        self.measurement_df = measurement_df

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
            self.pathway_dict = data['pathway_dict']

            # compound id -> formula, extracted from xml file
            self.entity_dict = data['entity_dict']

            # compound id -> [ mapids ], extracted from rdata file
            self.mapping_dict = data['mapping_dict']

            # map between pathway id to compound ids and formulas
            pathway_to_unique_ids_dict = defaultdict(set) # mapid -> [ formulas ]
            for entity_id, pathway_ids in self.mapping_dict.items():
                try:
                    unique_id = self.entity_dict[entity_id]['unique_id']
                    for pathway_id in pathway_ids:
                        pathway_to_unique_ids_dict[pathway_id].add(unique_id)
                except KeyError:
                    continue  # TODO: skip this compound since we can't find its name or formula. Should never happen!!
            self.pathway_to_unique_ids_dict = dict(pathway_to_unique_ids_dict)

            # a dataframe of peak id, originating database name, database id, formula
            dataset_pathways = []
            dataset_pathways_to_row_ids = defaultdict(list)
            for row_id, row in annotation_df.iterrows():
                entity_id = row['entity_id']
                try:
                    possible_pathways = self.mapping_dict[entity_id]
                    dataset_pathways.extend(possible_pathways)
                    for p in possible_pathways:
                        dataset_pathways_to_row_ids[p].append(row_id)
                except KeyError: # no information about compound -> pathway in our database
                    continue
            self.dataset_pathways = set(dataset_pathways)
            self.dataset_pathways_to_row_ids = dict(dataset_pathways_to_row_ids)
            self.dataset_unique_ids = set(annotation_df['unique_id'].values)

        # For use in the hypergeometric test - the number of unique formulas in kegg and in pathways
        # and the number of unique formulas in the ds and in pathways
        self.pathway_unique_ids_count = len(self.get_pathway_unique_ids())
        self.pathway_dataset_unique_ids_count = len(self.get_pathway_dataset_unique_ids())

    def get_pathway_unique_ids(self):
        """
        Returns the unique ids present in pathways
        :return: unique ids present in pathways
        """
        all_entity_ids = set(self.entity_dict.keys()) # all entity ids in database
        pathway_entity_ids = set(self.mapping_dict.keys()) # all entity ids found in pathways
        entity_ids_in_pathways = all_entity_ids.intersection(pathway_entity_ids)
        pathway_unique_ids = set([self.entity_dict[entity_id]['unique_id'] for entity_id in entity_ids_in_pathways])
        return pathway_unique_ids

    def get_pathway_dataset_unique_ids(self):
        """
        Returns the unique ids present in pathways in the dataset
        :return: unique ids present in pathways in the dataset
        """
        pathway_dataset_unique_ids = self.get_pathway_unique_ids().intersection(self.dataset_unique_ids)
        return pathway_dataset_unique_ids

    def get_pathway_unique_counts(self, pathway_ids):
        """
        Returns the number of unique ids associated to each pathway
        :param pathway_ids: pathway ids
        :return: a list of the number of unique ids associated to each pathway
        """
        return [len(self.pathway_to_unique_ids_dict[pathway_id]) for pathway_id in pathway_ids]

    def get_pathway_dataset_unique_counts(self, pathway_ids):
        """
        Gets the number of unique ids in the dataset for all the pathways
        :param pathway_ids: a list of pathway ids
        :return: a list of count of unique ids in the dataset for all the pathways
        """
        counts = []
        for pathway_id in pathway_ids:
            if pathway_id in self.dataset_pathways:  # get all the formulas of the pathway in the dataset
                unique_ids = self.pathway_to_unique_ids_dict[pathway_id].intersection(self.dataset_unique_ids)
                counts.append(len(unique_ids))
            else:
                counts.append(0)
        return counts

    def get_comparison_samples(self, comp):
        conditions = [comp['control'], comp['case']]
        condition_samples = list(map(lambda group_name: self.groups[group_name], conditions))
        return condition_samples
