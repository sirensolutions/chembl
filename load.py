import time
import json
from tqdm import tqdm
import os
import argparse
from elasticsearch import Elasticsearch, RequestsHttpConnection
from elasticsearch.helpers import streaming_bulk

tables = [
    'activities',
    'assays',
    'molecules',
    'papers',
    'target',
]

IMPORT_DIR = '/data/chembl/import_mols_only'

table2id = dict(
    activities='activity_id',
    assays='assay_id',
    molecules='molregno',
    papers='doc_id',
    target='tid'
)

class ChemblTableLoader():

    def __init__(self, es_url, es_user, es_pass, import_dir, table, limit=None, test=False):
        """

        :param es_url:
        :param es_user:
        :param es_pass:
        :param import_dir:
        :param table:
        :param limit:
        :param test:   add "test-" to start of index names so we don't overwrite the original
        """

        self.es = Elasticsearch(es_url,
                           verify_certs=False,
                           http_auth=(es_user, es_pass),
                           connection_class=RequestsHttpConnection)
        self.import_dir = import_dir
        self.table = table

        self.index_name = ('test-' if test else '') + 'chembl-%s' % table

        self.count = 0
        self.dump_file_name = os.path.join(self.import_dir, 'chembl-%s.json' % self.table)

    @staticmethod
    def get_numlines(dump_file_name):
        for i, line in enumerate(open(dump_file_name)):
            pass
        return i

    def data_iterator(self, id_field):
        num_lines = ChemblTableLoader.get_numlines(self.dump_file_name)
        print num_lines
        for i, line in tqdm(enumerate(open(self.dump_file_name)),
                            desc='loading %s in elasticsearch' % self.table,
                            total=num_lines):
            doc = json.loads(line)
            self.count += 1
            yield {
                '_index': self.index_name,
                '_type': 'document',
                '_id': doc[id_field],
                '_source': doc
            }

    def load_table_to_es(self):
        success, failed = 0, 0
        start_time = time.time()
        for ok, item in streaming_bulk(self.es, self.data_iterator(table2id[self.table]), raise_on_error=False, chunk_size=1000):
            if not ok:
                failed += 1
            else:
                success += 1
        print(
            'loading %s in es took %i seconds, %i success, %i failed ' % (self.table, time.time() - start_time, success, failed))

    def delete_index(self):
        print('deleting', self.index_name, self.es.indices.delete(index=self.index_name, ignore=404, timeout='300s'))

    def create_index(self):
        print('creating', self.index_name, self.es.indices.create(index=self.index_name, ignore=400, timeout='30s',
                                                        body=json.load(open('mappings/chembl-%s.json' % self.table))))



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Import chembl json into elasticsearch')
    parser.add_argument('-es', type=str, default='https://localhost:9220')
    parser.add_argument('-table', dest='tables', action='append', default=[], help='tables to load/update')
    parser.add_argument('-limit', type=int, default=None)
    parser.add_argument('-mappings', default='kibi')
    parser.add_argument('-username', required=False, default='admin')
    parser.add_argument('-password', required=False, default='password')
    parser.add_argument('-importdir', required=True)
    parser.add_argument('-test', action="store_true")

    args = parser.parse_args()

    for table in tables:
        table_loader = ChemblTableLoader(args.es, args.username, args.password, args.importdir, table, limit=args.limit, test=args.test)
        '''prepare indexes'''
        table_loader.delete_index()
        table_loader.create_index()
        '''load data'''
        table_loader.load_table_to_es()



