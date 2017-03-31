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


mappings are in the mappings directory but you can use a script we provide here.

 ./elasticsearchMapping.sh
Arguments:

	{Delete, Create}
 	index name
 	mapping file
 Example
	./elasticsearchMapping.sh Create chembl-activities mappings/chembl-activities.json
	./elasticsearchMapping.sh Delete chembl-activities

Create all the mappings that are in the mapping folder

4) Download logstash and make sure your JDBC plugin is installed.


Environment setup:
	Download and extract logstash:
		wget https://download.elastic.co/logstash/logstash/logstash-1.5.4.tar.gz
		tar zxf logstash-1.5.4.tar.gz
	Install logstash-input-jdbc
		cd logstash-1.5.4
		bin/plugin install logstash-input-jdbc
	Download jdbc driver
		wget https://jdbc.postgresql.org/download/postgresql-9.4-1204.jdbc42.jar
	Modify bin/logstash  by adding LS_HEAP_SIZE="2g" otherwise isues with heap size can occur

4) Get Kibi and start the Elasticsearch that comes with it. Feel free to start Kibi too.

Alternatively, use any Elasticsearch as long as it has the SIREn Join plugin http://siren.solutions

without the Siren Join plugins the relational features of the demo will not work.

5) Running logstash and create the Elasticsearch indexes

logstash-1.5.4/bin/logstash -f chembl-logstash/sql-target.jsn
logstash-1.5.4/bin/logstash -f chembl-logstash/sql-papers.jsn
logstash-1.5.4/bin/logstash -f chembl-logstash/sql-molecules.jsn
logstash-1.5.4/bin/logstash -f chembl-logstash/sql-assays.jsn
logstash-1.5.4/bin/logstash -f chembl-logstash/sql-activities.jsn

6) Now use Kibi (start if it you didnt start it before) and load the Dashboard configurations for this demo:

you can do so from the "settings/objects" menu in Kibi, pressing the "import" button and selecting the
kibi-configuration-objects.json file

if all went well you'll have all the dashboards and visualizations that are shown in the youtube video.

Notes:
* the SQL queries in the loader can be improved a LOT. e.g. there isnt much data about the molecules : just improve the SQL query that loads the molecules to use a "sister table" (linked 1:1 via molregno value ) that has much more infos about the molecules. Also paper have no authors or titles at the moment.
* the visualization of the molecules is pretty random :) the web services takes a parameter that is NOT the same as the molregno ID. It should be possible to get the equivalence somehow, if you know how to do it, please let us know :)
