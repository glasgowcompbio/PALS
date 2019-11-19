import gzip
import json
import os

PW_F_OFFSET = 1
MIN_REPLACE = 5000
NUM_RESAMPLES = 1000
PLAGE_WEIGHT = 5
HG_WEIGHT = 1

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
