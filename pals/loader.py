import os
import zipfile
from collections import defaultdict, Counter
from io import BytesIO

import pandas as pd
import requests
from loguru import logger
from tqdm import tqdm

from .common import DATABASE_PIMP_KEGG, load_json, DATA_DIR, GNPS_DOWNLOAD_CYTOSCAPE_DATA_VIEW, \
    GNPS_VIEW_ALL_MOTIFS_VIEW, \
    DATABASE_GNPS_MS2LDA, DATABASE_GNPS_MOLECULAR_FAMILY
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
    def __init__(self, database_name, gnps_url, metadata_df, comparisons, gnps_ms2lda_url=None, peak_table_df=None):
        self.database_name = database_name
        self.gnps_url = gnps_url
        self.metadata_df = metadata_df
        self.comparisons = comparisons
        self.int_df = None
        self.annotation_df = None
        self.gnps_ms2lda_url = gnps_ms2lda_url
        self.peak_table_df = peak_table_df

        if self.database_name == DATABASE_GNPS_MS2LDA:
            assert self.gnps_ms2lda_url is not None

    def load_data(self):
        if self.peak_table_df is not None:  # load measurements from a peak table
            logger.info('Processing peak table')
            logger.debug(self.peak_table_df)

            # drop the first (m/z) and second (RT) columns to get the measurement df
            cols = [0, 1]
            measurement_df = self.peak_table_df.drop(self.peak_table_df.columns[cols], axis=1)
            measurement_df.index.rename('peak_id', inplace=True)
            measurement_df.index = measurement_df.index.astype('str')

            # FIXME: really shouldn't be called this
            clustering_df = self.peak_table_df[self.peak_table_df.columns[cols]]

            # create annotation dataframe
            annotation_df = pd.DataFrame(index=measurement_df.index)
            annotation_df.index.rename('peak_id', inplace=True)
            annotation_df['entity_id'] = measurement_df.index
            annotation_df['entity_id'] = annotation_df['entity_id'].astype(str)
            annotation_df.index = annotation_df.index.astype('str')

        else:  # load measurements from GNPS
            logger.info('Retrieving clustering and quantification information from GNPS')
            logger.debug(self.gnps_url)
            results = self._download_gnps(self.gnps_url, GNPS_DOWNLOAD_CYTOSCAPE_DATA_VIEW)
            assert results is not None
            quantification_df = results['quantification_df']
            clustering_df = results['clustering_df']
            filtered_clustering_df = clustering_df[
                clustering_df['componentindex'] != -1]  # drop all the singleton components

            # keep only columns containing 'Peak area', and remove 'Peak area' from column names
            measurement_df = quantification_df.filter(regex='Peak area')
            measurement_df.columns = measurement_df.columns.str.rstrip('Peak area')
            measurement_df.index.rename('peak_id', inplace=True)
            measurement_df.index = measurement_df.index.astype('str')

            # create annotation dataframe
            annotation_df = pd.DataFrame(index=filtered_clustering_df.index)
            annotation_df.index.rename('peak_id', inplace=True)
            annotation_df['entity_id'] = filtered_clustering_df.index
            annotation_df['entity_id'] = annotation_df['entity_id'].astype(str)
            annotation_df.index = annotation_df.index.astype('str')

        # filter dataframes
        # assume metadata_df has two columns: 'sample' and 'group'
        # remove rows with sample id that can't be found in the columns of int_df
        metadata_df = self.metadata_df[self.metadata_df['sample'].isin(measurement_df.columns.values)]
        # keep only columns in int_df that have group information
        measurement_df = measurement_df[metadata_df['sample']]

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

        # Turn grouping information into PALS database object
        # If it's a standard FBMN-GNPS result, then use the clustering as the groups
        # otherwise if it is GNPS-MS2LDA result, then download the MS2LDA results from GNPS and use motifs as groups
        if self.database_name == DATABASE_GNPS_MOLECULAR_FAMILY:
            filtered_clustering_df = filtered_clustering_df.rename(columns={
                'precursor mass': 'mass',
                'RTConsensus': 'RT'
            })
            database = self._molfam_to_database(filtered_clustering_df, extra_data)

        elif self.database_name == DATABASE_GNPS_MS2LDA:
            logger.info('Retrieving motif information from GNPS')
            logger.debug(self.gnps_ms2lda_url)

            results = self._download_gnps(self.gnps_ms2lda_url, GNPS_VIEW_ALL_MOTIFS_VIEW)
            motif_df = results['motif_df']

            # select some useful columns to display later
            # need to include all peaks, so we select the columns from clustering_df (instead of filtered_clustering_df)
            try:
                peak_info_df = clustering_df[['parent mass', 'LibraryID', 'GNPSLinkout_Network', 'number of spectra',
                                              'RTConsensus', 'precursor mass', 'SumPeakIntensity', 'componentindex']]
                peak_info_df = peak_info_df.rename(columns={
                    'precursor mass': 'mass',
                    'RTConsensus': 'RT'
                })
            except KeyError:
                peak_info_df = clustering_df[['mass', 'RT']]
            database = self._motif_to_database(peak_info_df, motif_df, extra_data)

        return database

    def _molfam_to_database(self, clustering_df, extra_data):
        """
        Creates a user-defined database from GNPS Molecular Family clustering
        :param clustering_df: a dataframe of GNPS clustering information
        :param extra_data: additional information to include in the database
        :return: a Database object from GNPS Molecular Family clustering
        """

        # Create 'pathway' dictionary. In this case, 'pathway' is a GNPS molecular family
        pathway_dict = {}
        for comp in clustering_df['componentindex'].values:
            key = str(comp)
            pathway_dict[key] = {'display_name': 'Molecular Family #%d' % comp}

        # Create entity dictionary. An 'entity' is a MS1 peak (GNPS consensus cluster)
        entity_dict = self._get_entity_dict(clustering_df)

        # Create mapping dictionary that maps entities to pathways
        mapping_dict = {}
        for peak_id in entity_dict:
            component_index = str(entity_dict[peak_id]['componentindex'])
            mapping_dict[peak_id] = [component_index]

        # put everything together in a Database object
        database = Database(self.database_name, pathway_dict, entity_dict, mapping_dict, extra_data=extra_data)
        return database

    def _motif_to_database(self, peak_info_df, motif_df, extra_data):
        """
        Creates a user-defined database from GNPS-MS2LDA results
        :param peak_info_df: a dataframe of additional information for peaks
        :param motif_df: a dataframe of LDA analysis from GNPS-MS2LDA
        :param extra_data: additional information to include in the database
        :return: a Database object from GNPS-MS2LDA results
        """
        # find singleton motifs
        c = Counter()
        for idx, row in motif_df.iterrows():
            motif = row['motif']
            c[motif] += 1

        motifs = motif_df['motif'].unique()
        singletons = [motif for motif in motifs if c[motif] == 1]

        # Create 'pathway' dictionary. In this case, 'pathway' is a GNPS-MS2LDA motif
        pathway_dict = {}
        motifdb_urls = {}
        for idx, row in motif_df.iterrows():
            key = row['motif']
            if key in singletons:
                continue

            motifdb_url = row['motifdb_url']
            motifdb_annotation = row['motifdb_annotation']

            # Try to cast motifdb_annotation to float. If success, then it contains NaN, which we can ignore
            # otherwise add motifdb_annotation to the display name
            try:
                float(motifdb_annotation)  # will throw ValueError if this contains an annotation string
                display_name = key
            except ValueError:
                display_name = '%s [%s]' % (key, motifdb_annotation)

            pathway_dict[key] = {
                'display_name': '%s' % display_name,
                'motifdb_url': motifdb_url,
                'motifdb_annotation': motifdb_annotation
            }
            motifdb_urls[display_name] = motifdb_url

        # Create entity dictionary. An 'entity' is a MS1 peak (GNPS consensus cluster)
        entity_dict = self._get_entity_dict(peak_info_df)

        # Create mapping dictionary that maps entities to pathways
        mapping_dict = defaultdict(list)
        for idx, row in motif_df.iterrows():
            peak_id = str(row['scan'])
            motif = row['motif']
            if motif in singletons:
                continue
            mapping_dict[peak_id].append(motif)
        mapping_dict = dict(mapping_dict)

        # put everything together in a Database object
        extra_data['motifdb_urls'] = motifdb_urls
        database = Database(self.database_name, pathway_dict, entity_dict, mapping_dict, extra_data=extra_data)
        return database

    def _get_entity_dict(self, peak_info_df):
        # First turn the peak info dataframe to dictionary, with peak id as the key
        peak_info_df.index = peak_info_df.index.astype('str')
        temp = peak_info_df.to_dict(orient='index')

        # Extract entity information from temp
        # temp contains a lot of stuff we don't want, so copy selected values to entity_dict
        entity_dict = {}
        for peak_id in temp:
            entity_dict[peak_id] = {}
            entity_dict[peak_id]['unique_id'] = peak_id
            entity_dict[peak_id]['mass'] = temp[peak_id]['mass']
            entity_dict[peak_id]['RT'] = temp[peak_id]['RT']
            try:
                entity_dict[peak_id]['display_name'] = temp[peak_id]['parent mass']
                entity_dict[peak_id]['LibraryID'] = temp[peak_id]['LibraryID']
                entity_dict[peak_id]['GNPSLinkout_Network'] = temp[peak_id]['GNPSLinkout_Network']
                entity_dict[peak_id]['number of spectra'] = temp[peak_id]['number of spectra']
                entity_dict[peak_id]['SumPeakIntensity'] = temp[peak_id]['SumPeakIntensity']
                entity_dict[peak_id]['componentindex'] = temp[peak_id]['componentindex']
            except KeyError:
                pass
        return entity_dict

    def _download_gnps(self, gnps_url, view):
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
            'view': view
        }
        api_endpoint = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult'
        r = requests.post(url=api_endpoint, data=data, stream=True)

        # extract clustering and quantification tables
        # https://stackoverflow.com/questions/37573483/progress-bar-while-download-file-over-http-with-requests
        total_size = int(r.headers.get('content-length', 0))
        block_size = 1024
        results = None
        with BytesIO() as f, tqdm(total=total_size, unit='iB', unit_scale=True) as t:
            for data in r.iter_content(block_size):
                t.update(len(data))
                f.write(data)

            if view == GNPS_DOWNLOAD_CYTOSCAPE_DATA_VIEW:
                clustering_df, quantification_df = self._parse_gnps_molfam(f)
                results = {
                    'clustering_df': clustering_df,
                    'quantification_df': quantification_df
                }
            elif view == GNPS_VIEW_ALL_MOTIFS_VIEW:
                motif_df = self._parse_ms2lda_motifs(f)
                results = {
                    'motif_df': motif_df
                }
        return results

    def _parse_gnps_molfam(self, input_stream):
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

    def _parse_ms2lda_motifs(self, input_stream):
        motif_df = None
        with zipfile.ZipFile(input_stream) as z:
            for filename in z.namelist():
                if 'view_all_motifs' in filename:
                    logger.debug('Found motif table: %s' % filename)
                    motif_df = pd.read_csv(z.open(filename), sep='\t')

        return motif_df
