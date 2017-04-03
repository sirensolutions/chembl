#!/usr/bin/env bash
mkdir -p export
rm -f export/*.json
#kibi
elasticdump --input=http://localhost:9201/.kibi --type=mapping --output=export/kibi-mappings.json
elasticdump --input=http://localhost:9201/.kibi --type=data --output=export/kibi-data.json
#activities
elasticdump --input=http://localhost:9201/chembl-activities --type=analyzer --output=export/chembl-activities-analyzer.json
elasticdump --input=http://localhost:9201/chembl-activities --type=mapping --output=export/chembl-activities-mappings.json
#assays
elasticdump --input=http://localhost:9201/chembl-assays --type=analyzer --output=export/chembl-assays-analyzer.json
elasticdump --input=http://localhost:9201/chembl-assays --type=mapping --output=export/chembl-assays-mappings.json
#papers
elasticdump --input=http://localhost:9201/chembl-papers --type=analyzer --output=export/chembl-papers-analyzer.json
elasticdump --input=http://localhost:9201/chembl-papers --type=mapping --output=export/chembl-papers-mappings.json
#target
elasticdump --input=http://localhost:9201/chembl-target --type=analyzer --output=export/chembl-target-analyzer.json
elasticdump --input=http://localhost:9201/chembl-target --type=mapping --output=export/chembl-target-mappings.json
#molecules
elasticdump --input=http://localhost:9201/chembl-molecules --type=analyzer --output=export/chembl-molecules-analyzer.json
elasticdump --input=http://localhost:9201/chembl-molecules --type=mapping --output=export/chembl-molecules-mappings.json

cp import/*.json export

#compress
tar -cvzf export/chembk-kibi-data.tar.gz export/*.json