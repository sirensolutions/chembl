SIREN INVESTIGATE LIFE SCIENCE DEMO - CHEMBL
======================

The code in this project can be used to create ChEMBL indexes as found in the life demo.


1) Start a local fingerprint API (Docker required).
    ```
    cd fingerprint_api
    docker build -t fp_api:latest .
    # docker run -d --name fp_api fp_api:latest
    docker run --name ls-api -p 8009:8080  -e PORT=8080 --rm fp_api
    ```
    
    OR, to allow the script to compute the fingerprints locally
    
    ```apt-get install -y python-pip python-dev build-essential python-rdkit librdkit1 rdkit-data```

2) Create a venv to run the scripts
   
   ```python -m virtualenv  --system-site-packages  -p /usr/bin/python venv```
   ```source venv/bin/activate```

3) Download and extract the latest _sqlite.tar.gz file, e.g. 

    ```ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latestchembl_24_1_sqlite.tar.gz```

4) Run extract_chembl.py 

    ```python extract_chembl.py -importdir chembl/import -chemblsqlite chembl_24/chembl_24_sqlite/chembl_24.db```

    this will:
   * Query the database and generate the input json files in the 'importdir' directory
   * Digest all the compounds into compatible fingerprints

5) Run load.py to load the json into elastic

    ```python load.py -importdir chembl/import -es http://localhost:9201  -api https://localhost:8009``` 

    using the urls your elasticsearch instance and (optionally)the fingerprint API are exposed to and wait for it to complete
   
   * Load the json in elasticsearch with the proper mapping
   * Create the necessary index-pattern in Investigate
