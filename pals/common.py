import gzip
import json
import logging
import os
import pathlib
import pickle
import sys

from loguru import logger

PW_F_OFFSET = 1
MIN_REPLACE = 5000
NUM_RESAMPLES = 1000
PLAGE_WEIGHT = 5
HG_WEIGHT = 1
SIGNIFICANT_THRESHOLD = 0.05

DATABASE_PIMP_KEGG = 'PiMP_KEGG'
DATABASE_REACTOME_KEGG = 'COMPOUND'
DATABASE_REACTOME_CHEBI = 'ChEBI'
DATABASE_REACTOME_UNIPROT = 'UniProt'
DATABASE_REACTOME_ENSEMBL = 'ENSEMBL'

REACTOME_SPECIES_ARABIDOPSIS_THALIANA = 'Arabidopsis thaliana'
REACTOME_SPECIES_BOS_TAURUS = 'Bos taurus'
REACTOME_SPECIES_CAENORHABDITIS_ELEGANS = 'Caenorhabditis elegans'
REACTOME_SPECIES_CANIS_LUPUS_FAMILIARIS = 'Canis lupus familiaris'
REACTOME_SPECIES_DANIO_RERIO = 'Danio rerio'
REACTOME_SPECIES_DICTYOSTELIUM_DISCOIDEUM = 'Dictyostelium discoideum'
REACTOME_SPECIES_DROSOPHILA_MELANOGASTER = 'Drosophila melanogaster'
REACTOME_SPECIES_GALLUS_GALLUS = 'Gallus gallus'
REACTOME_SPECIES_HOMO_SAPIENS = 'Homo sapiens'
REACTOME_SPECIES_MUS_MUSCULUS = 'Mus musculus'
REACTOME_SPECIES_ORYZA_SATIVA = 'Oryza sativa'
REACTOME_SPECIES_RATTUS_NORVEGICUS = 'Rattus norvegicus'
REACTOME_SPECIES_SACCHAROMYCES_CEREVISIAE = 'Saccharomyces cerevisiae'
REACTOME_SPECIES_SUS_SCROFA = 'Sus scrofa'

PATHWAY_ANALYSIS_PALS = 'PALS'
PATHWAY_ANALYSIS_ORA = 'ORA'

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def load_json(json_file, compressed=False):
    if compressed:
        with gzip.GzipFile(json_file, 'r') as f:
            data = json.loads(f.read().decode('utf-8'))
    else:
        with open(json_file, 'r') as f:
            data = json.load(f)
    return data


def save_json(data, json_file, compressed=False):
    if compressed:
        with gzip.GzipFile(json_file, 'w') as f:
            f.write(json.dumps(data).encode('utf-8'))
    else:
        with open(json_file, 'w') as f:
            json.dump(data, f)


def create_if_not_exist(out_dir):
    if not os.path.exists(out_dir) and len(out_dir) > 0:
        logger.debug('Created %s' % out_dir)
        pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)


def save_obj(obj, filename):
    """
    Save object to file
    :param obj: the object to save
    :param filename: the output file
    :return: None
    """
    out_dir = os.path.dirname(filename)
    create_if_not_exist(out_dir)
    logger.debug('Saving %s to %s' % (type(obj), filename))
    with gzip.GzipFile(filename, 'w') as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_obj(filename):
    """
    Load saved object from file
    :param filename: The file to load
    :return: the loaded object
    """
    try:
        with gzip.GzipFile(filename, 'rb') as f:
            return pickle.load(f)
    except OSError:
        logger.warning('Old, invalid or missing pickle in %s. Please regenerate this file.' % filename)
        return None


def set_log_level_warning():
    logger.remove()
    logger.add(sys.stderr, level=logging.WARNING)


def set_log_level_info():
    logger.remove()
    logger.add(sys.stderr, level=logging.INFO)


def set_log_level_debug():
    logger.remove()
    logger.add(sys.stderr, level=logging.DEBUG)


def is_comparison_used(comp, selected_case, selected_control):
    # skip this comparison if not the same as what the user has specified
    comp_case = comp['case']
    comp_control = comp['control']
    if selected_case is not None and selected_control is not None:
        if selected_case != comp_case or selected_control != comp_control:
            return False
    return True
