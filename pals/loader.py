import os
from io import BytesIO

import pandas as pd
import requests
import zipfile
from loguru import logger
from tqdm import tqdm

from .common import DATABASE_PIMP_KEGG, load_json, DATA_DIR
from .reactome import get_pathway_dict, get_compound_mapping_dict, load_entity_dict, get_protein_entity_dict, \
    get_protein_mapping_dict, get_gene_entity_dict, get_gene_mapping_dict


class Database(object):
    def __init__(self, database_name, pathway_dict, entity_dict, mapping_dict, extra_data=None):
        self.database_name = database_name
        self.pathway_dict = pathway_dict
        self.entity_dict = entity_dict
        self.mapping_dict = mapping_dict
        self.extra_data = extra_data

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


class GNPSLoader(Loader):
    def __init__(self, database_name, gnps_url, metadata_df, comparisons):
        self.database_name = database_name
        self.gnps_url = gnps_url
        self.metadata_df = metadata_df
        self.comparisons = comparisons
        self.int_df = None
        self.annotation_df = None

    def load_data(self):
        # load clustering and quantification info from the zip file
        clustering_df, quantification_df = self._download_gnps(self.gnps_url)
        assert clustering_df is not None and quantification_df is not None

        # drop all the singleton components
        clustering_df = clustering_df[clustering_df['componentindex'] != -1]

        # keep only columns containing 'Peak area', and remove 'Peak area' from column names
        measurement_df = quantification_df.filter(regex='Peak area')
        measurement_df.columns = measurement_df.columns.str.rstrip('Peak area')
        measurement_df.index.rename('peak_id', inplace=True)
        measurement_df.index = measurement_df.index.astype('str')

        # assume metadata_df has two columns: 'sample' and 'group'
        # remove rows with sample id that can't be found in the columns of int_df
        metadata_df = self.metadata_df[self.metadata_df['sample'].isin(measurement_df.columns.values)]

        # keep only columns in int_df that have group information
        measurement_df = measurement_df[metadata_df['sample']]

        # create annotation dataframe
        annotation_df = pd.DataFrame(index=clustering_df.index)
        annotation_df.index.rename('peak_id', inplace=True)
        annotation_df['entity_id'] = clustering_df.index
        annotation_df['entity_id'] = annotation_df['entity_id'].astype(str)
        annotation_df.index = annotation_df.index.astype('str')

        # create experimental design dictionary
        groups = {}
        for k, v in metadata_df.groupby('group'):
            groups[k] = v['sample'].values.tolist()

        experimental_design = {
            'comparisons': self.comparisons,
            'groups': groups
        }

        # combine all above into the extra_data dictionary for a Database
        extra_data = {
            'measurement_df': measurement_df,
            'annotation_df': annotation_df,
            'experimental_design': experimental_design
        }

        # generate a database for GNPS
        database = self._to_database(clustering_df, extra_data)
        return database

    def _to_database(self, clustering_df, extra_data):
        """
        Creates a user-defined database for GNPS
        :param clustering_df: a dataframe of clustering information
        :param extra_data: additional information to include in the database
        :return: a Database object for GNPS
        """

        # First create 'pathway' dictionary. In this case, 'pathway' is a GNPS molecular family / component
        pathway_dict = {}
        for comp in clustering_df['componentindex'].values:
            key = str(comp)
            pathway_dict[key] = {'display_name': 'Molecular Family #%d' % comp}

        # Create entity dictionary. In this case, an 'entity' is a MS1 peak (GNPS consensus cluster)
        # turn df to dictionary, with cluster index as the key
        clustering_df.index = clustering_df.index.astype('str')
        entity_dict = clustering_df.to_dict(orient='index')
        for peak_id in entity_dict:
            entity_dict[peak_id]['unique_id'] = peak_id
            entity_dict[peak_id]['display_name'] = entity_dict[peak_id]['parent mass']

        # Create mapping dictionary that maps entities to pathways
        mapping_dict = {}
        for peak_id in entity_dict:
            component_index = str(entity_dict[peak_id]['componentindex'])
            mapping_dict[peak_id] = [component_index]

        # put everything together in a Database object
        database = Database(self.database_name, pathway_dict, entity_dict, mapping_dict, extra_data=extra_data)
        return database

    def _download_gnps(self, gnps_url):
        """
        Downloads the zipped cytoscape data from GNPS and extract clustering and quantification dataframes from it.
        :param gnps_url: the url to the GNPS experiment, e.g.
        https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=0a8432b5891a48d7ad8459ba4a89969f
        :return: clustering and quantification dataframes from that GNPS result
        """
        # extract task id from the url
        tokens = gnps_url.split('task=')
        task = tokens[1]
        logger.info('Found GNPS task %s' % task)

        # send a post request to GNPS
        data = {
            'task': task,
            'view': 'download_cytoscape_data'
        }
        api_endpoint = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult'
        r = requests.post(url=api_endpoint, data=data, stream=True)

        # extract clustering and quantification tables
        # https://stackoverflow.com/questions/37573483/progress-bar-while-download-file-over-http-with-requests
        clustering_df = None
        quantification_df = None
        total_size = int(r.headers.get('content-length', 0))
        block_size = 1024

        t = tqdm(total=total_size, unit='iB', unit_scale=True)
        with BytesIO() as f:
            for data in r.iter_content(block_size):
                t.update(len(data))
                f.write(data)
            clustering_df, quantification_df = self._parse_gnps(f)
        t.close()

        return clustering_df, quantification_df

    def _parse_gnps(self, input_stream):
        """
        Parses a zipped GNPS input stream, and extract clustering and quantification tables
        :param input_stream: a zipped input of GNPS results
        :return: clustering and quantification tables from the zip file
        """
        clustering_df = None
        quantification_df = None
        with zipfile.ZipFile(input_stream) as z:
            for filename in z.namelist():
                # read clustering information
                if filename.startswith('clusterinfo_summary'):
                    logger.debug('Found cluster info: %s' % filename)
                    clustering_df = pd.read_csv(z.open(filename), sep='\t', index_col='cluster index')

                # read quantification information
                if filename.startswith('quantification_table'):
                    logger.debug('Found quantification table: %s' % filename)
                    quantification_df = pd.read_csv(z.open(filename), sep=',').set_index('row ID')

        return clustering_df, quantification_df
