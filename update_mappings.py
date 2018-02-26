import argparse
import elasticsearch


'''import input'''
parser = argparse.ArgumentParser(description='Import mappings into elasticsearch')
parser.add_argument('-es', type=str, default='https://localhost:9220')
parser.add_argument('-table', dest='tables', action='append',
                    default=[],
                    help='tables to update the mappings')

args = parser.parse_args()

tables = [
    'activities',
    'assays',
    'molecules',
    'papers',
    'target',
]

for table in tables:
    '''prepare indexes'''
    index_name = 'chembl-%s' % table
    print('deleting', index_name, es.indices.delete(index=index_name, ignore=404, timeout='300s'))
    print('creating', index_name, es.indices.create(index=index_name, ignore=400, timeout='30s',
                                                    body=json.load(open('mappings/%s.json' % index_name))))
