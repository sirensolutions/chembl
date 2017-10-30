import argparse
import json
import os
import sqlite3
import time
import warnings
from subprocess import call

import numpy as np
import requests
from elasticsearch import Elasticsearch, RequestsHttpConnection
from elasticsearch.helpers import parallel_bulk, streaming_bulk
from tqdm import tqdm

warnings.filterwarnings("ignore")

'''import input'''
parser = argparse.ArgumentParser(description='Import chembl data into elasticsearch')
parser.add_argument('-es', type=str, default='https://localhost:9222')
args = parser.parse_args()

'''SETUP'''
CHEMBL_DB_VERSION = 'chembl_23'
ES_URL = args.es
ES_AUTH = ('admin', 'password')
CHEMBL_SQLITE_DB_DIR = CHEMBL_DB_VERSION + '_sqlite'
CHEMBL_SQLITE_DB = os.path.join(CHEMBL_SQLITE_DB_DIR, CHEMBL_DB_VERSION + '.db')
CHEMBL_DB_DUMP_FILE = CHEMBL_SQLITE_DB_DIR + '.tar.gz'
CHEMBL_SQLITE_URL = 'http://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/%s/%s' % (
    CHEMBL_DB_VERSION, CHEMBL_DB_DUMP_FILE)
CONCAT_SEPARATOR = '|'
IMPORT_DIR = 'import'
if not os.path.exists(IMPORT_DIR):
    os.mkdir(IMPORT_DIR)

s = requests.Session()
'''download database file'''
# if not os.path.exists(CHEMBL_DB_DUMP_FILE) or \
#         not os.path.exists(CHEMBL_SQLITE_DB):
#     print CHEMBL_DB_DUMP_FILE, CHEMBL_SQLITE_DB
#     # r = requests.get(CHEMBL_SQLITE_URL, stream=True)
#     # total_size = int(r.headers.get('content-length', 0));
#     #
#     # with open(CHEMBL_SQLITE_DB_DIR+'.tar.gz', 'wb') as f:
#     #     for data in tqdm(r.iter_content(32*1024),
#     #                      total=total_size,
#     #                      unit='B',
#     #                      unit_scale=True,
#     #                      desc='Download database dump'):
#     #         f.write(data)
#     '''download file'''
#     call(["curl", '--output', CHEMBL_DB_DUMP_FILE, '-O', CHEMBL_SQLITE_URL])
#     '''uncompress file'''
#     call(["tar", "zxvf", CHEMBL_DB_DUMP_FILE])


def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d


tables = ['activities',
          'assays',
          'molecules',
          'papers',
          'target',
          ]

queries = dict(activities='''SELECT
  activities.molregno,
  activities.assay_id,
  activities.activity_id,
  activities.doc_id,
  activities.standard_relation,
  activities.standard_value,
  activities.standard_units,
  activities.standard_flag,
  activities.standard_type,
  activities.activity_comment,
  ((lower(activity_comment) LIKE '%not active%') OR (lower(activity_comment) LIKE '%inactive%'))  AS inactive,
  activities.data_validity_comment
FROM
  activities''',
               assays='''SELECT
  assays.assay_id,
  assays.assay_type,
  assay_type.assay_type,
  assay_type.assay_desc,
  assays.description,
  assays.assay_category,
  assays.tid,
  assays.confidence_score,
  assays.chembl_id,
  assays.doc_id,
  assays.relationship_type,
  assays.assay_test_type,
  assays.assay_organism,
  assays.assay_strain,
  assays.assay_tissue,
  assays.assay_subcellular_fraction,
  relationship_type.relationship_desc
FROM
  assays,
  assay_type,
  relationship_type
WHERE
  assays.assay_type = assay_type.assay_type AND
  relationship_type.relationship_type = assays.relationship_type;
  ''',
               molecules='''SELECT
  group_concat(molecule_synonyms.synonyms, '{0}') AS synonyms,
  molecule_dictionary.pref_name,
  molecule_dictionary.molregno,
  molecule_dictionary.chembl_id,
  molecule_dictionary.therapeutic_flag,
  molecule_dictionary.molecule_type,
  molecule_dictionary.chirality,
  molecule_dictionary.inorganic_flag,
  molecule_dictionary.polymer_flag,
  molecule_dictionary.indication_class,
  molecule_dictionary.structure_type,
  molecule_dictionary.usan_year,
  molecule_dictionary.availability_type,
  compound_properties.*,
  biotherapeutics.description,
  biotherapeutics.helm_notation,
  drug_indication.max_phase_for_ind,
  group_concat(drug_indication.efo_id, '{0}') AS efo_id,
  group_concat(drug_indication.efo_term, '{0}') AS efo_term,
  group_concat(drug_indication.mesh_id, '{0}') AS mesh_id,
  group_concat(drug_indication.mesh_heading, '{0}') AS mesh_heading,
  compound_structures.canonical_smiles,
  cr.compound_name,
  cr.compound_doc_id,
  cr.compound_source_description,
  cr.src_short_name
FROM
  molecule_dictionary
  LEFT JOIN (SELECT
      compound_records.molregno,
      group_concat(compound_records.compound_name, '{0}') AS compound_name,
      group_concat(compound_records.doc_id, '{0}') AS compound_doc_id,
      group_concat(source.src_description, '{0}') AS compound_source_description,
      group_concat(source.src_short_name, '{0}') AS src_short_name
    FROM
      compound_records
      LEFT JOIN source ON compound_records.src_id = source.src_id
    GROUP BY compound_records.molregno) as cr
      ON molecule_dictionary.molregno = cr.molregno
  LEFT JOIN molecule_synonyms ON molecule_synonyms.molregno = molecule_dictionary.molregno
  LEFT JOIN compound_structures ON molecule_dictionary.molregno = compound_structures.molregno
  LEFT JOIN compound_properties ON molecule_dictionary.molregno = compound_properties.molregno
  LEFT JOIN biotherapeutics ON molecule_dictionary.molregno = biotherapeutics.molregno
  LEFT JOIN drug_indication On molecule_dictionary.molregno = drug_indication.molregno
GROUP BY molecule_dictionary.molregno
  '''.format(CONCAT_SEPARATOR),
               papers='''SELECT
  docs.doc_id,
  docs.journal,
  docs.year,
  docs.volume,
  docs.issue,
  docs.first_page,
  docs.pubmed_id,
  docs.last_page,
  docs.doi,
  docs.chembl_id,
  docs.title,
  docs.authors,
  docs.abstract,
  docs.doc_type
FROM
  docs
  ''',
               target='''SELECT
  target_dictionary.tid,
  target_dictionary.pref_name,
  target_dictionary.organism,
  target_dictionary.chembl_id,
  target_dictionary.target_type,
  group_concat(component_synonyms.component_synonym,'%s') AS synonyms,
  component_sequences.accession

FROM
  target_dictionary
  LEFT JOIN target_components ON target_components.tid = target_dictionary.tid
  LEFT JOIN component_synonyms ON component_synonyms.component_id = target_components.component_id
  LEFT JOIN component_sequences ON component_sequences.component_id = target_components.component_id
GROUP BY target_dictionary.tid;
  ''' % CONCAT_SEPARATOR)

table2id = dict(activities='activity_id',
                assays='assay_id',
                molecules='molregno',
                papers='doc_id',
                target='tid'
                )


def encode_vector(v):
    e = []
    for i, x in enumerate(v):
        e.append(str(i) if x > 0 else 'z' + str(i))
    return e


'''Export data to json files'''
try:
    db = sqlite3.connect(CHEMBL_SQLITE_DB)
    db.row_factory = dict_factory  # sqlite3.Row
    cursor = db.cursor()
except:
    print 'no chembl database available'

extracted_counts = dict()
for table in tables:
    start_time = time.time()
    i = 0
    dump_file_name = os.path.join(IMPORT_DIR, 'chembl-%s.json' % table)
    if not os.path.exists(dump_file_name):
        cursor.execute(queries[table])
        print 'Extracting data for table %s' % table
        with open(dump_file_name, 'w') as f:
            for i, row in tqdm(enumerate(cursor)):
                for k, v, in row.items():
                    try:
                        if CONCAT_SEPARATOR in v:
                            row[k] = list(set(v.split(CONCAT_SEPARATOR)))
                    except TypeError:
                        pass
                if table == 'molecules' and row['canonical_smiles']:
                    fingerprint_r = s.get('https://localhost:8009/binaryfingerprint',
                                          params={'smiles':row['canonical_smiles']},
                                          verify=False)
                    if fingerprint_r.ok:
                        fingerprint = fingerprint_r.json()
                        # fingerprint_int = fingerprint.astype(np.int8).tolist()
                        # row['fingerprint'] = fingerprint_int
                        # row['fingerprint_b'] = (fingerprint>0).astype(np.int8).tolist()
                        # row['fingerprint_nz'] = ' '.join([str(j) for j in fingerprint.nonzero()[0].tolist()])
                        row['fingerprint_all'] = ' '.join(encode_vector(map(int,list(fingerprint))))
                        # row['fingerprint_minhash'] = row['fingerprint_nz']
                    else:
                        print i,'error finger print for smiles: '+ row['canonical_smiles']

                f.write(json.dumps(row) + '\n')
        print('exporting table %s took %i seconds, %i rows' % (table, time.time() - start_time, i))
    else:
        for i, line in enumerate(open(dump_file_name)):
            pass
    extracted_counts[table] = i

'''Import json files in elasticsearch'''
es = Elasticsearch(ES_URL,
                   use_ssl=True,
                   verify_certs=False,
                   http_auth=ES_AUTH,
                   connection_class=RequestsHttpConnection)


def data_iterator(table, id_field):
    for i, line in tqdm(enumerate(open(os.path.join(IMPORT_DIR, 'chembl-%s.json' % table))),
                        desc='loading %s in elasticsearch' % table,
                        total=extracted_counts[table]):
        doc = json.loads(line)
        yield {
            '_index': 'chembl-%s' % table,
            '_type': 'document',
            '_id': doc[id_field],
            '_source': doc
        }


def load_table_to_es(table):
    success, failed = 0, 0
    start_time = time.time()
    for ok, item in streaming_bulk(es,
                                  data_iterator(table,
                                                table2id[table]),
                                  raise_on_error=False,
                                  chunk_size=1000):
        if not ok:
            failed += 1
        else:
            success += 1
    print(
        'loading %s in es took %i seconds, %i success, %i failed ' % (table, time.time() - start_time, success, failed))


for table in tables:
    '''prepare indexes'''
    index_name = 'chembl-%s' % table
    print('deleting', index_name, es.indices.delete(index=index_name, ignore=404, timeout='300s'))
    print('creating', index_name, es.indices.create(index=index_name, ignore=400, timeout='30s',
                                                    body=json.load(open('mappings/%s.json' % index_name))))

    '''load data'''
    load_table_to_es(table)

# '''load empty kibi configuration'''
# for table in tables:
#     index_name = 'chembl-%s' % table
#     es.index(index='.kibi',
#              doc_type='index-pattern',
#              id=index_name,
#              body={"title": index_name})
# es.index(index='.kibi',
#          doc_type='config',
#          id="4.6.4",
#          body={"defaultIndex": "chembl-activities"})


# '''load kibi preconfired index'''
# index_name = '.kibi'
# print('deleting', index_name, es.indices.delete(index=index_name, ignore=404, timeout='300s'))
# print('creating', index_name, es.indices.create(index=index_name, ignore=400, timeout='30s',
#                                                 body=json.load(open('kibi/kibi-mappings.json'))['.kibi']))
#
# '''load objects'''
# success, failed = 0, 0
# for ok, item in streaming_bulk(es,
#                                (json.loads(i) for i in open('kibi/kibi-data.json').readlines()),
#                                raise_on_error=False,
#                                chunk_size=1000):
#     if not ok:
#         failed += 1
#     else:
#         success += 1
#
# print('loaded %i objects in .kibi index. %i failed' % (success, failed))
