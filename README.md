SIREN INVESTIGATE LIFE SCIENCE DEMO
======================

The code in this repo can reproduce the latest Siren Investigate ChEMBL demo.
These are the steps to reproduce the demo.

1) Start up a Investigate instance. E.g. using docker:
    ```docker run  -d -p 5606:5606 -p 9201:9220 --net=chembl --name chembl-kibi sirensolutions/siren-platform:5.4.3-2```

2) [Optional] Start a local fingerprint API. Docker required.
    ```
    cd fingerprint_api
    docker build -t fp_api:latest .
    # docker run -d --name fp_api fp_api:latest
    docker run --name ls-api -p 8009:8080  -e PORT=8080 --rm fp_api
    ```

3) ```pip install -r requirements.txt```

4) run ```python import.py -es http://localhost:9201  -api https://localhost:8009``` using the urls your elasticsearch
instance and the fingerprint API are exposed to and wait for it to complete
   this will:
   * Download and extract the chembl dumps in SQLite format
   * Query the database and generate the input json file in the `import` directory
   * Digest all the compounds into compatible fingerprints
   * Load the json in elasticsearch with the proper mapping
   * Create the necessary index-pattern in Investigate

5) explore the data at your Investigate instance: http://localhost:5606

```
$ cd ~/src/repos/kibi-internal
$ grep 9220 config/investigate.yml
elasticsearch.url: "http://localhost:9220"
$ npm start
```

6) To persist the state of Investigate you can run `dump_kibi.sh`. This requires [elasticdump](https://www.npmjs.com/package/elasticdump) to be installed.
   To save the data and the kibi configuartion in a compressed file you can run `dump_all.sh`.


6 Alternative ways to persist the data:

```
$ # To bypass unauthorised connections
$ export NODE_TLS_REJECT_UNAUTHORIZED=0
$ elasticdump --input=https://<username>:<password>@localhost:9220/.kibi --output=kibi/kibi-mappings.json —type=mapping
$ elasticdump --input=https://<username>:<password>@localhost:9220/.kibiaccess --output=kibi/kibiacces-mappings.json —type=mapping
$ elasticdump --input=https://<username>:<password>@localhost:9220/.kibiaccess --output=kibi/kibiacces-data.json —type=data
```

Alternatively you can dump all the files with `bin/investigate backup`.

```
$ # To bypass unauthorised connections
$ export NODE_TLS_REJECT_UNAUTHORIZED=0
$ bin/investigate backup --dev --backup-dir investigate
```


NOTES
-----
* This demo requires python >2.7.10 in order for the TLS1.2 to work with the searchguard plugin installed in kibi >5.

* The fingerprint API can be quickly deployed on google app engine using the prepared `app.yaml` file

* To include ensembl ids in the `chembl-target` index, download the mappings (from Ensembl's BioMart or UniProt), format the file in CSV having the uniprot IDs in the first column and the Ensembl Ids in the second and run `include_ensembl_ids.py` script.
