import argparse
import warnings
from collections import defaultdict
from tqdm import tqdm

from elasticsearch import Elasticsearch, RequestsHttpConnection
from elasticsearch.helpers import scan


warnings.filterwarnings('ignore')

'''import input'''
parser = argparse.ArgumentParser(description='Import Uniprot - Ensembl mappings into chembl-papers')
parser.add_argument('-es', type=str, default='https://localhost:9220')
parser.add_argument('-mappings', type=str, default='./mart_export_uniprot_ensembl_human.tsv')
parser.add_argument('-username', type=str, required=False, default='admin')
parser.add_argument('-password', type=str, required=False, default='password')
args = parser.parse_args()

def get_number_of_records(es):
    # query = {
    #     'query': {}
    # }
    #
    r = es.count(index='chembl-target')
    return r['count']


if __name__ == '__main__':
    es = Elasticsearch(args.es,
                       use_ssl=True,
                       verify_certs=False,
                       http_auth=(args.username, args.password),
                       connection_class=RequestsHttpConnection)

    mapping_dict = defaultdict(list)
    '''Read mapping file and import into es'''
    with open(args.mappings, 'r') as mappings:
        for cnt, line in enumerate(mappings):
            if cnt == 0:
                continue
            else:
                parts = line.rstrip().split('\t')
                if len(parts) > 1:
                    mapping_dict[parts[1]].append(parts[0])

        c = get_number_of_records(es)
        with tqdm(
            desc='loading ensembl ids in targets',
            unit=' docs',
            unit_scale=True,
            total=c
        ) as pbar:
            for item in scan(client=es,
                               query={
                                   'query': {}
                               },
                               scroll='5m',
                               raise_on_error=True,
                               preserve_order=False,
                               size=100,
                               index='chembl-target'):
                pbar.update(1)

                # print(item)
                uniprotId = item['_source']['accession']
                doc_id = item['_id']
                doc_type = item['_type']
                # print(uniprotId)
                # print(mapping_dict[uniprotId])
                es.update(index='chembl-target',
                          id=doc_id,
                          doc_type=doc_type,
                          body={
                              'doc': {
                                  'ensembl_ids': mapping_dict[uniprotId]
                              }
                          })
