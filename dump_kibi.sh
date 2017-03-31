#!/usr/bin/env bash
mkdir -p export
rm -f export/kibi*.json
#kibi
elasticdump --input=http://localhost:9201/.kibi --type=mapping --output=export/kibi-mappings.json
elasticdump --input=http://localhost:9201/.kibi --type=data --output=export/kibi-data.json
