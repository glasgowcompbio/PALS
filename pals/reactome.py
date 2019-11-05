import os
import re
from collections import defaultdict

from neo4j import GraphDatabase, basic_auth

from common import DATA_DIR, DATABASE_PIMP_KEGG, load_json, save_json, DATABASE_REACTOME_KEGG

NEO4J_SERVER = 'bolt://localhost:7687'
NEO4J_USER = 'neo4j'
NEO4J_PASSWORD = 'neo4j'
driver = GraphDatabase.driver(NEO4J_SERVER,
                              auth=basic_auth(NEO4J_USER, NEO4J_PASSWORD))


def get_neo4j_session():
    session = None
    try:
        session = driver.session()
    except Exception:
        raise
    return session


def rchop(thestring, ending):
    if thestring.endswith(ending):
        return thestring[:-len(ending)]
    return thestring


def get_species_list():
    results = []
    try:
        session = get_neo4j_session()
        query = """
        MATCH (n:Species) RETURN n.displayName AS name order by name        
        """
        query_res = session.run(query)
        # print(query)
        for record in query_res:
            results.append(record['name'])
    finally:
        if session is not None: session.close()
    return results


def get_pathway_dict(species, metabolic_pathway_only=True, leaf=True):
    results = {}
    try:
        session = get_neo4j_session()

        # initial match clause in the query
        query = """
            MATCH (tp:TopLevelPathway)-[:hasEvent*]->(p:Pathway)-[:hasEvent*]->(rle:ReactionLikeEvent)
            WHERE
                tp.speciesName = {species} AND
        """

        if leaf:  # retrieve only the leaf nodes in the pathway hierarchy
            query += " (p)-[:hasEvent]->(rle) AND "

        if metabolic_pathway_only:  # only retrieves metabolic pathways
            query += " tp.displayName = 'Metabolism' AND "

        # remove last AND
        query = rchop(query.strip(), 'AND')

        # add return clause
        query += """
            RETURN DISTINCT
                p.speciesName AS species_name,            
                p.displayName AS pathway_name,
                p.stId AS pathway_id                       
        """

        params = {
            'species': species
        }
        query_res = session.run(query, params)
        # print(query)

        for record in query_res:
            pathway_name = record['pathway_name']
            pathway_id = record['pathway_id']
            results[pathway_id] = {'display_name': pathway_name}
    finally:
        if session is not None: session.close()
    return results


def get_compound_mapping_dict(species, database_name, metabolic_pathway_only=True, leaf=True):
    results = defaultdict(list)
    try:
        session = get_neo4j_session()

        # initial match clause in the query
        query = """
        MATCH (tp:TopLevelPathway)-[:hasEvent*]->
              (p:Pathway)-[:hasEvent*]->(rle:ReactionLikeEvent),
              (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent
              |hasMember|hasCandidate*]->(pe:PhysicalEntity),
              (pe:PhysicalEntity)-[:crossReference|:referenceEntity]->(do:DatabaseObject)
        WHERE
              tp.speciesName = {species} AND
              do.databaseName = {database_name} AND
        """

        if leaf:  # retrieve only the leaf nodes in the pathway hierarchy
            query += " (p)-[:hasEvent]->(rle) AND "

        if metabolic_pathway_only:  # only retrieves metabolic pathways
            query += " tp.displayName = 'Metabolism' AND "

        # remove last AND
        query = rchop(query.strip(), 'AND')

        # add return clause
        query += """
            RETURN DISTINCT
                p.stId AS pathway_id,
                do.identifier AS entity_id
        """

        params = {
            'species': species,
            'database_name': database_name
        }
        query_res = session.run(query, params)
        # print(query)

        i = 0
        for record in query_res:
            pathway_id = record['pathway_id']
            entity_id = record['entity_id']
            results[entity_id].append(pathway_id)
    finally:
        if session is not None: session.close()
    return dict(results)


def write_database(pathway_dict, entity_dict, mapping_dict, json_file):
    if len(mapping_dict) > 0:
        data = {
            'pathway_dict': pathway_dict,
            'entity_dict': entity_dict,
            'mapping_dict': mapping_dict
        }
        save_json(data, json_file, compressed=True)


def load_entity_dict(database_name):
    json_file = os.path.join(DATA_DIR, '%s.json.zip' % database_name)
    return load_json(json_file, compressed=True)


def parse_chebi_entity_dict(owl_file):
    """
    Parses ChEBI compound dict from the Ontology file downloaded from https://www.ebi.ac.uk/chebi/downloadsForward.do
    :param owl_file: the ChEBI ontology file
    :return: a ChEBI entity dictionary
    """
    chebi_id = None
    display_name = None
    formula = None
    chebi_entity_dict = {}

    with open(owl_file, encoding='utf-8') as f:
        for line in f:
            if 'owl:Class' in line and 'rdf:about' in line:
                found = line.strip()
                res = re.search('CHEBI_(.*)"', found)
                chebi_id = res.group(1)
            if 'chebi:formula' in line:
                found = line.strip()
                res = re.search('<chebi:formula.*>(.*)<\/chebi:formula>', found)
                formula = res.group(1)
            if 'rdfs:label' in line:
                found = line.strip()
                res = re.search('<rdfs:label.*>(.*)<\/rdfs:label>', found)
                display_name = res.group(1)
            if '</owl:Class>' in line:
                if chebi_id is not None and display_name is not None and formula is not None:
                    chebi_entity_dict[chebi_id] = {
                        'display_name': display_name,
                        'unique_id': formula
                    }
                    chebi_id = None
                    display_name = None
                    formula = None
    return chebi_entity_dict
