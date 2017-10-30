KIBI ChEMBL DEMO II
===================

The code in this repo can reproduce the second Kibi ChEMBL demo.
This are the steps to reproduce the demo.

1) Start up a kibi instance. E.g. using docker:
    ```docker run  -d -p 5606:5606 -p 9201:9220 --net=chembl --name chembl-kibi sirensolutions/kibi-community-standard:4.6.4```

2) ```pip install -r requirements.txt```

3) run ```python import.py -es http://localhost:9201``` using the port your elasticsearch instance is exposed to and wait for it to complete
   this will:
   * Download and extract the chembl dumps in SQLite format
   * Query the database and generate the input json file in the `import` directory
   * Load the json in elasticsearch with the proper mapping
   * Create the necessary index-pattern in kibi

4) explore the data at your kibi instance: http://localhost:5606

5) To persist the state of kibi you can run `dump_kibi.sh`. This requires [elasticdump](https://www.npmjs.com/package/elasticdump) to be installed.
   To save the data and the kibi configuartion in a compressed file you can run `dump_all.sh`.

NOTES
-----
This demo requires python >2.7.10 in order for the TLS1.2 to work with the searchguard plugin installed in kibi >5