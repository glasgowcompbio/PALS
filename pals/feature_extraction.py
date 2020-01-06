import os
from collections import defaultdict

import numpy as np
import pandas as pd
from loguru import logger
from sklearn import preprocessing

from .common import DATABASE_PIMP_KEGG, load_json, DATA_DIR, DATABASE_REACTOME_KEGG, DATABASE_REACTOME_CHEBI, \
    DATABASE_REACTOME_UNIPROT, DATABASE_REACTOME_ENSEMBL, MIN_REPLACE
from .reactome import get_pathway_dict, get_compound_mapping_dict, load_entity_dict, get_protein_entity_dict, \
    get_protein_mapping_dict, get_gene_entity_dict, get_gene_mapping_dict


class DataSource(object):

    def __init__(self, measurement_df, annotation_df, experimental_design, database_name,
                 reactome_species=None, reactome_metabolic_pathway_only=True, reactome_query=False, database=None,
                 min_replace=MIN_REPLACE):
        """
        Creates a data source for PALS analysis
        :param measurement_df: a dataframe of peak intensities, where index = row id and columns = sample_name
        :param annotation_df: a dataframe where index = row id. There are two required columns:
               - entity_id is the annotation assigned to row id, e.g. kegg compound
               - unique_id is used to uniquely indicate which entity it is, e.g. formula of the kegg compound
        :param experimental_design: a dictionary specifying the experimental design
        :param database_name: the database name (with .json extension) used to load database in the data folder
        """
        self.measurement_df = measurement_df
        self.annotation_df = annotation_df
        self.experimental_design = experimental_design
        self.database_name = database_name
        self.reactome_species = reactome_species
        self.reactome_metabolic_pathway_only = reactome_metabolic_pathway_only
        self.reactome_query = reactome_query
        self.min_replace = min_replace

        self.groups = dict(self.experimental_design['groups'].items())
        self.comparisons = self.experimental_design['comparisons']

        if database is None:
            logger.debug('Using %s as database' % database_name)
            assert database_name is not None
            self.database = self.get_database(database_name, reactome_metabolic_pathway_only, reactome_query,
                                              reactome_species)
        else:
            logger.debug('Using user-provided database')
            assert database_name is None
            self.database = database
        self.pathway_dict = self.database['pathway_dict']  # mapid -> pathway name
        self.entity_dict = self.database['entity_dict']  # compound id -> formula
        self.mapping_dict = self.database['mapping_dict']  # compound id -> [ mapids ]

        # map between pathway id to compound ids and formulas
        logger.debug('Mapping pathway to unique ids')
        pathway_to_unique_ids_dict = defaultdict(set)  # mapid -> [ formulas ]
        for entity_id, pathway_ids in self.mapping_dict.items():
            try:
                # if unique_id is present, then use it
                # useful when we want to represent chemicals as formulae
                unique_id = self._get_unique_id(entity_id)
                for pathway_id in pathway_ids:
                    pathway_to_unique_ids_dict[pathway_id].add(unique_id)
            except KeyError:
                continue  # TODO: skip this compound since we can't find its name or formula. Should never happen!!
        self.pathway_to_unique_ids_dict = dict(pathway_to_unique_ids_dict)

        # a dataframe of peak id, originating database name, database id, formula
        logger.debug('Creating dataset to pathway mapping')
        dataset_pathways = []
        dataset_pathways_to_row_ids = defaultdict(list)
        dataset_unique_ids = []
        dataset_row_id_to_unique_ids = defaultdict(set)
        for row_id, row in annotation_df.iterrows():
            entity_id = row['entity_id']
            try:
                # convert entity id to unique id if possible
                unique_id = self._get_unique_id(entity_id)
                dataset_row_id_to_unique_ids[row_id].add(unique_id)
                dataset_unique_ids.append(unique_id)
                # get the mapping between dataset pathway to row ids
                possible_pathways = self.mapping_dict[entity_id]
                dataset_pathways.extend(possible_pathways)
                for p in possible_pathways:
                    dataset_pathways_to_row_ids[p].append(row_id)
            except KeyError:
                # no unique id found for this entity id, OR
                # no information about compound -> pathway in our database
                continue
        self.dataset_pathways = set(dataset_pathways)
        self.dataset_pathways_to_row_ids = dict(dataset_pathways_to_row_ids)
        self.dataset_unique_ids = set(dataset_unique_ids)
        self.dataset_row_id_to_unique_ids = dataset_row_id_to_unique_ids

        # For use in the hypergeometric test - the number of unique formulas in kegg and in pathways
        # and the number of unique formulas in the ds and in pathways
        logger.debug('Computing unique id counts')
        self.pathway_unique_ids_count = len(self._get_pathway_unique_ids())
        self.pathway_dataset_unique_ids_count = len(self._get_pathway_dataset_unique_ids())

    def get_database(self, database_name, reactome_metabolic_pathway_only, reactome_query, reactome_species):
        # load compound and pathway database information from file
        if database_name == DATABASE_PIMP_KEGG:  # PiMP exported pathways for KEGG
            json_file = os.path.abspath(os.path.join(DATA_DIR, '%s.json.zip' % DATABASE_PIMP_KEGG))
            logger.debug('Loading %s' % json_file)
            data = load_json(json_file, compressed=True)

        elif database_name == DATABASE_REACTOME_KEGG or database_name == DATABASE_REACTOME_CHEBI:  # must be using reactome
            if reactome_query:  # fetch reactome data from neo4j
                logger.debug('Retrieving data for %s from Reactome %s metabolic_pathway_only=%s' %
                             (reactome_species, database_name, reactome_metabolic_pathway_only))
                pathway_dict = get_pathway_dict(reactome_species, reactome_metabolic_pathway_only)
                mapping_dict = get_compound_mapping_dict(reactome_species, database_name,
                                                         reactome_metabolic_pathway_only)
                entity_dict = load_entity_dict(database_name)
                data = {
                    'pathway_dict': pathway_dict,
                    'entity_dict': entity_dict,
                    'mapping_dict': mapping_dict
                }
            else:
                # we didn't dump the data for all pathways. Only for the metabolic pathways only this can be used.
                if not reactome_metabolic_pathway_only:
                    raise ValueError(
                        'Pathway information is not available. Please use live reactome query with --connect_to_reactome_server.')
                metabolic_pathway_dir = 'metabolic_pathways' if reactome_metabolic_pathway_only else 'all_pathways'
                json_file = os.path.join(DATA_DIR, 'reactome', metabolic_pathway_dir, database_name,
                                         '%s.json.zip' % reactome_species)
                logger.debug('Loading %s' % json_file)
                data = load_json(json_file, compressed=True)

        elif database_name == DATABASE_REACTOME_UNIPROT:
            pathway_dict = get_pathway_dict(reactome_species, reactome_metabolic_pathway_only)
            entity_dict = get_protein_entity_dict(reactome_species, database_name)
            mapping_dict = get_protein_mapping_dict(reactome_species, database_name)
            data = {
                'pathway_dict': pathway_dict,
                'entity_dict': entity_dict,
                'mapping_dict': mapping_dict
            }

        elif database_name == DATABASE_REACTOME_ENSEMBL:
            pathway_dict = get_pathway_dict(reactome_species, reactome_metabolic_pathway_only)
            entity_dict = get_gene_entity_dict(reactome_species, database_name)
            mapping_dict = get_gene_mapping_dict(reactome_species, database_name)
            data = {
                'pathway_dict': pathway_dict,
                'entity_dict': entity_dict,
                'mapping_dict': mapping_dict
            }
        return data

    def get_measurements(self):
        return self.measurement_df.copy()

    def get_annotations(self):
        return self.annotation_df.copy()

    def get_experimental_design(self):
        return self.experimental_design.copy()

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
        """
        Given a comparison, e.g. {'case': 'beer1', 'control': 'beer2', 'name': 'beer1/beer2'}
        returns the samples from the control and the case (in that order)
        :param comp: a dictionary that describes a comparison
        :return: a list containing two lists. The first is the control samples, and the next
        is the case samples. Example:
        [['Beer_2_full3.mzXML', 'Beer_2_full1.mzXML', 'Beer_2_full2.mzXML'],
         ['Beer_1_full2.mzXML', 'Beer_1_full1.mzXML', 'Beer_1_full3.mzXML']]
        """
        conditions = [comp['control'], comp['case']]
        condition_samples = list(map(lambda group_name: self.groups[group_name], conditions))
        return condition_samples

    def resample(self, n_sample, case=None, control=None, axis=0):
        """
        Resamples columns for performance evaluation
        :param case: the case label
        :param control: the control label
        :param n_sample: the number of columns to sample
        :return: a new DataSource object containing the new resampled data
        """
        measurement_df = self.get_measurements()
        annotation_df = self.get_annotations()
        if axis == 0:
            sampled_measurement_df = measurement_df.sample(n_sample, axis=0, replace=False)
            sampled_annotation_df = annotation_df[annotation_df.index.isin(sampled_measurement_df.index)]
            experimental_design = self.get_experimental_design()
            new_ds = DataSource(sampled_measurement_df, sampled_annotation_df, experimental_design, self.database_name,
                                self.reactome_species, self.reactome_metabolic_pathway_only, self.reactome_query)

        elif axis == 1:
            assert case is not None, 'Case is required'
            assert control is not None, 'Control is required'
            samples_case = self.groups[case]
            samples_control = self.groups[control]

            # sample n_samples columns without replacement
            intensities_case = measurement_df[samples_case].sample(n_sample, axis=1, replace=False)
            intensities_control = measurement_df[samples_control].sample(n_sample, axis=1, replace=False)

            # combined the sampled case and control dataframes together
            combined_df = pd.concat([intensities_case, intensities_control], axis=1)

            # sample the entire dataframe again to randomise the column order
            shuffled_df = combined_df.sample(frac=1, axis=1)

            # create a filtered experimental design
            selected_samples = set(shuffled_df.columns.values.tolist())
            experimental_design = {
                'comparisons': [{
                    'case': case,
                    'control': control,
                    'name': '%s/%s' % (case, control)
                }],
                'groups': {
                    case: list(set(samples_case).intersection(selected_samples)),
                    control: list(set(samples_control).intersection((selected_samples)))
                }
            }

            # return a new DataSource
            new_ds = DataSource(shuffled_df, annotation_df, experimental_design, self.database_name,
                                self.reactome_species, self.reactome_metabolic_pathway_only, self.reactome_query)
        return new_ds

    def standardize_intensity_df(self):
        """
        Standardize measurement dataframe by filling in missing values and standardizing across samples
        :param measurement_df: Dataframe of measured intensites (raw)
        :return: DF with zero intensities replaced and the values standardized
        """
        # Change the 0.00 intensities in the matrix to useable values
        measurement_df = self.change_zero_peak_ints()

        # standardize the data across the samples (zero mean and unit variance))
        logger.debug("Scaling the data across the sample: zero mean and unit variance")
        scaled_data = np.log(np.array(measurement_df))
        scaled_data = preprocessing.scale(scaled_data, axis=1)

        # Put the scaled data back into df for further use
        sample_names = measurement_df.columns
        measurement_df[sample_names] = scaled_data
        return measurement_df

    def change_zero_peak_ints(self):
        """
        A method to change a 'zero' entries in a dataframe.
        If all intensities in a (factor) group are zero, a min value is set.
        If there are > 1 and < number in group zero intensities, then the average of the non_zeros entries is calculated
        and used. Assuming the PiMP mzXML file names are unique
        :param measurement_df: A dataframe of peak intensities with peak ids (rows) and samples (columns)
        :return: No return, modifies peak_int_df.
        """
        # Get the min_intensity value set for the analysis
        logger.debug("Setting the zero intensity values in the dataframe")
        measurement_df = self.get_measurements()

        # Replace 0.0 with NaN for easier operations ahead
        measurement_df[measurement_df == 0.0] = None
        for group_name, samples in self.groups.items():
            # If all zero in group then replace with minimum
            measurement_df.loc[measurement_df.loc[:, samples].isnull().all(axis=1), samples] = self.min_replace

            # Replace any other zeros with mean of group
            subset_df = measurement_df.loc[:, samples]
            measurement_df.loc[:, samples] = subset_df.mask(subset_df.isnull(), subset_df.mean(axis=1), axis=0)
        return measurement_df

    def _calculate_coverage_df(self, mapids):
        """
        Calculate the Formula coverage for a dataset.
        :param mapids: The Mapids for the pathways to be used
        :return: A dataframe containing the number of unique formulae for a pathway, along with those
        annotated, identified and the total unique Fomulae in each pathway.
        """
        logger.debug("Calculating dataset formula coverage")

        # num_formula: Stores the number of unqique kegg formulae for a pathway
        num_formula = self.get_pathway_unique_counts(mapids)
        num_totalF = self.get_pathway_dataset_unique_counts(mapids)

        # unq_pw_f: unique formula expected for a pathway
        # tot_ds_F: unique formula for a pathway in a dataset
        data = {'unq_pw_F': num_formula,
                'tot_ds_F': num_totalF, }
        for_df = pd.DataFrame(data)
        for_df.index = mapids

        # Calculate the coverage of the formula found in the ds vs formulae in the pathway
        for_df['F_coverage'] = (((for_df['tot_ds_F']) / for_df['unq_pw_F']) * 100).round(2)
        return for_df

    def _get_unique_id(self, entity_id):
        """
        Returns unique id of an entity if 'unique_id' is present in the entity dictionary.
        Useful when we want to represent kegg chemical id as a compound
        :param entity_id: the unique id
        :return: either the unique id if available, or the original entity id otherwise
        """
        if 'unique_id' in self.entity_dict[entity_id]:
            unique_id = self.entity_dict[entity_id]['unique_id']
        else:
            unique_id = entity_id
        return unique_id

    def _get_pathway_unique_ids(self):
        """
        Returns the unique ids present in pathways
        :return: unique ids present in pathways
        """
        all_entity_ids = set(self.entity_dict.keys())  # all entity ids in database
        pathway_entity_ids = set(self.mapping_dict.keys())  # all entity ids found in pathways
        entity_ids_in_pathways = all_entity_ids.intersection(pathway_entity_ids)
        pathway_unique_ids = set([self._get_unique_id(entity_id) for entity_id in entity_ids_in_pathways])
        return pathway_unique_ids

    def _get_pathway_dataset_unique_ids(self):
        """
        Returns the unique ids present in pathways in the dataset
        :return: unique ids present in pathways in the dataset
        """
        pathway_dataset_unique_ids = self._get_pathway_unique_ids().intersection(self.dataset_unique_ids)
        return pathway_dataset_unique_ids
