from __future__ import print_function

from elasticsearch import Elasticsearch, RequestsHttpConnection
from flask import Flask, jsonify, request
from flask_cors import CORS
import chemfp
from rdkit import Chem
from rdkit.Chem import AllChem


app = Flask(__name__)
CORS(app)




def build_query(field,like, size=10):
    query = {
        "query": {
            "more_like_this": {
                "fields": [],
                "like": "",
                "min_term_freq": 1,
                "max_query_terms": 50,
                "minimum_should_match": "85%",
                "boost_terms": 2
            }
        },
        "_source": ["chembl_id", "compound_name", "canonical_smiles"],
        "size": size
    }
    query['query']['more_like_this']['fields'] = [field]
    query['query']['more_like_this']['like'] = like
    return query

def encode_vector(v):
    e = []
    for i, x in enumerate(v):
        e.append(str(i) if x > 0 else 'z' + str(i))
    return e

def get_fingerprint_from_smiles(smiles):
    m = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048, useFeatures=0, useChirality=0, useBondTypes=1)
    return fp.ToBitString()

def get_query_from_smiles(smiles):

    encoded_fingerprint = encode_vector(map(int, get_fingerprint_from_smiles(smiles)))
    # return build_query(field='fingerprint_all',
    #                            like=encoded_fingerprint)
    return ' '.join(encoded_fingerprint)



@app.route('/fingerprint')
def encoded_rdkit_morgan_fingerprint():
    smiles = request.args.get('smiles')
    return jsonify(get_query_from_smiles(smiles))


@app.route('/binaryfingerprint')
def rdkit_morgan_fingerprint():
    smiles = request.args.get('smiles')
    return jsonify(get_fingerprint_from_smiles(smiles))




if __name__ == '__main__':
    # to run locally  in https with self signed certificates:
    # context = ('12597520-api.life.siren.io.cert', '12597520-api.life.siren.io.key')
    # app.run(host='0.0.0.0', port=8008,
    #         ssl_context= context)
    app.run(host='0.0.0.0', port=8080)
