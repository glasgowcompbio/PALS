import os

from loguru import logger

from .common import DATABASE_PIMP_KEGG, load_json, DATA_DIR
from .reactome import get_pathway_dict, get_compound_mapping_dict, load_entity_dict, get_protein_entity_dict, \
    get_protein_mapping_dict, get_gene_entity_dict, get_gene_mapping_dict


class Database(object):
    def __init__(self, database_name, pathway_dict, entity_dict, mapping_dict):
        self.database_name = database_name
        self.pathway_dict = pathway_dict
        self.entity_dict = entity_dict
        self.mapping_dict = mapping_dict

    def __repr__(self):
        return self.database_name


class Loader(object):
    def load_data(self):
        raise NotImplementedError()


class PiMP_KEGG_Loader(Loader):
    def __init__(self, database_name):
        self.database_name = database_name

    def load_data(self):
        json_file = os.path.abspath(os.path.join(DATA_DIR, '%s.json.zip' % DATABASE_PIMP_KEGG))
        logger.debug('Loading %s' % json_file)
        data = load_json(json_file, compressed=True)
        database = Database(self.database_name, data['pathway_dict'], data['entity_dict'], data['mapping_dict'])
        return database


class CompoundOnlineLoader(Loader):
    def __init__(self, database_name, reactome_species, mp_only):
        self.database_name = database_name
        self.reactome_species = reactome_species
        self.mp_only = mp_only

    def load_data(self):
        logger.debug('Retrieving data for %s from Reactome %s metabolic_pathway_only=%s' %
                     (self.reactome_species, self.database_name, self.mp_only))
        pathway_dict = get_pathway_dict(self.reactome_species,
                                        metabolic_pathway_only=self.mp_only)
        mapping_dict = get_compound_mapping_dict(self.reactome_species, self.database_name,
                                                 metabolic_pathway_only=self.mp_only)
        entity_dict = load_entity_dict(self.database_name)
        data = {
            'pathway_dict': pathway_dict,
            'entity_dict': entity_dict,
            'mapping_dict': mapping_dict
        }
        database = Database(self.database_name, data['pathway_dict'], data['entity_dict'], data['mapping_dict'])
        return database


class CompoundOfflineLoader(Loader):
    def __init__(self, database_name, reactome_species, mp_only):
        self.database_name = database_name
        self.reactome_species = reactome_species
        self.mp_only = mp_only

    def load_data(self):
        if not self.mp_only:
            raise ValueError(
                'Pathway information is not available. Please use live reactome query with --connect_to_reactome_server.')
        metabolic_pathway_dir = 'metabolic_pathways' if self.mp_only else 'all_pathways'
        json_file = os.path.join(DATA_DIR, 'reactome', metabolic_pathway_dir, self.database_name,
                                 '%s.json.zip' % self.reactome_species)
        logger.debug('Loading %s' % json_file)
        data = load_json(json_file, compressed=True)
        database = Database(self.database_name, data['pathway_dict'], data['entity_dict'], data['mapping_dict'])
        return database


class UniProtLoader(Loader):
    def __init__(self, database_name, reactome_species, mp_only):
        self.database_name = database_name
        self.reactome_species = reactome_species
        self.mp_only = mp_only

    def load_data(self):
        pathway_dict = get_pathway_dict(self.reactome_species, metabolic_pathway_only=self.mp_only)
        entity_dict = get_protein_entity_dict(self.reactome_species, self.database_name)
        mapping_dict = get_protein_mapping_dict(self.reactome_species, self.database_name,
                                                metabolic_pathway_only=self.mp_only)
        data = {
            'pathway_dict': pathway_dict,
            'entity_dict': entity_dict,
            'mapping_dict': mapping_dict
        }
        database = Database(self.database_name, data['pathway_dict'], data['entity_dict'], data['mapping_dict'])
        return database


class EnsemblLoader(Loader):
    def __init__(self, database_name, reactome_species, mp_only):
        self.database_name = database_name
        self.reactome_species = reactome_species
        self.mp_only = mp_only

    def load_data(self):
        pathway_dict = get_pathway_dict(self.reactome_species, metabolic_pathway_only=self.mp_only)
        entity_dict = get_gene_entity_dict(self.reactome_species, self.database_name)
        mapping_dict = get_gene_mapping_dict(self.reactome_species, self.database_name,
                                             metabolic_pathway_only=self.mp_only)
        data = {
            'pathway_dict': pathway_dict,
            'entity_dict': entity_dict,
            'mapping_dict': mapping_dict
        }
        database = Database(self.database_name, data['pathway_dict'], data['entity_dict'], data['mapping_dict'])
        return database
