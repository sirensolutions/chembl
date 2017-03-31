KIBI ChEMBL DEMO II
===================

This demo works on the sqlite dump of Chembl.
A python script take queries the database to generate json objects and push them to elasticsearch.
It also generate the required index patterns so you can start exploring the data with kibi.

Steps are

1) Start up a kibi instance. E.g. using docker:
    ```docker run  -d -p 5606:5606 -p 9201:9220 --net=chembl --name chembl-kibi sirensolutions/kibi-community-standard:4.6.4```

2) ```pip install -r requirements.txt```

3) run ```python import.py -es http://localhost:9201``` using the port your elasticsearch isntance is exposed to and wait for it to complete

4) explore the data at your kibi instance: `http://localhost:5606'
