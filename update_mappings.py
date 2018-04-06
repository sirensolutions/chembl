import argparse
import json
import warnings
from elasticsearch import Elasticsearch, RequestsHttpConnection
from elasticsearch.helpers import parallel_bulk, streaming_bulk

warnings.filterwarnings("ignore")


'''import input'''
parser = argparse.ArgumentParser(description='Import mappings into elasticsearch')
parser.add_argument('-es', type=str, default='https://localhost:9220')
parser.add_argument('-mappings', default='kibi')
parser.add_argument('-username', required=False, default='admin')
parser.add_argument('-password', required=False, default='password')
args = parser.parse_args()

ES_URL = args.es
ES_AUTH = (args.username, args.password)

es = Elasticsearch(ES_URL,
                   use_ssl=True,
                   verify_certs=False,
                   http_auth=ES_AUTH,
                   connection_class=RequestsHttpConnection)

index_name = '.siren'
print('deleting', index_name, es.indices.delete(index=index_name, ignore=404, timeout='300s'))
print('deleting', '.kibi', es.indices.delete(index='.kibi', ignore=404, timeout='300s'))
print('deleting', '.investigate', es.indices.delete(index='.investigate', ignore=404, timeout='300s'))
print('creating', index_name, es.indices.create(index=index_name, ignore=400, timeout='30s',
                                                body=json.load(open('%s/mapping-.kibi.json' % (args.mappings)))[
                                                    '.siren']))

index_name = '.sirenaccess'
print('deleting', index_name, es.indices.delete(index=index_name, ignore=404, timeout='300s'))
print('deleting', '.kibiaccess', es.indices.delete(index='.kibiaccess', ignore=404, timeout='300s'))
print('deleting', '.investigateaccess', es.indices.delete(index='.investigateaccess', ignore=404, timeout='300s'))
print('creating', index_name, es.indices.create(index=index_name, ignore=400, timeout='30s',
                                                body=json.load(open('%s/mapping-.kibiaccess.json' % (args.mappings)))[
                                                    '.sirenaccess']))

success, failed = 0, 0
for ok, item in streaming_bulk(es,
                               (json.loads(i) for i in open('%s/data-.kibi.json' % (args.mappings)).readlines()),
                               raise_on_error=False,
                               chunk_size=1000):
    if not ok:
        failed += 1
    else:
        success += 1

print('loaded %i objects in .siren index. %i failed' % (success, failed))

success, failed = 0, 0
for ok, item in streaming_bulk(es,
                               (json.loads(i) for i in open('%s/data-.kibiaccess.json' % (args.mappings)).readlines()),
                               raise_on_error=False,
                               chunk_size=1000):
    if not ok:
        failed += 1
    else:
        success += 1

print('loaded %i objects in .sirenaccess index. %i failed' % (success, failed))
