Metab Import
============

## Prerequisites

Setup a Python Virtual Environment, then install the required dependencies.

    $ python3 -m venv venv
    $ source venv/bin/activate
    $ pip install -r requirements.txt


## Preparing your config file

To run the importer, you will need a config file containing two sections: one related to the Postgres database, and one related to accessing VIVO. The Postgres section includes credentials for accessing the database, the open port you used to set up the SSH tunnel, and the name of the database. The VIVO section includes credentials for accessing the VIVO query API, as well as the query endpoint and the namespace used for VIVO URIs. There is an example config file with the required fields. Download a CSV export from the [Metabolomics Software Tools](https://docs.google.com/spreadsheets/d/1a096jlzbAwTxUdvTtJB2bwfkBRdDcMKrkywS8YZf344/edit?usp=sharing) and note the file name in the config for 'csv_tools'.


## Connect to the database

The first step to using the importer is setting up an SSH tunnel to stage.vivo.metabolomics.info to access the Postgres database with the information from Metabolomics Workbench. You will need to set up the tunnel in its own terminal tab.

    $ ssh -N -L 5432:localhost:5432 stage.vivo.metabolomics.info


## Running the Pre-fill script

Once the tunnel is prepared, you can run the pre-fill script.

    $ python metab_prefill.py $CONFIG_PATH

This pre-fills the supplemental tables with necessary information like people
and organizations.


## Running the Importer

Next, run:

    $ python metab_import.py $CONFIG_PATH

This will produce up to four files: projects.rdf, studies.rdf, datasets.rdf, and people.rdf. These files contain the triples for each respective class. After the files are printed, the importer will delete the entire http://vitro.mannlib.cornell.edu/default/vitro-kb-2 graph on your VIVO instance. **Everything in this graph will be removed**. If you have manually added triples that are not a part of the import process, they *will* be lost.

If you wish to print the files without deleting the database or automatically uploading the importer output, you can use the `-d` or `--dry-run` flag.

    $ python metab_import.py -d $CONFIG_PATH

## Running the Admin Page

To start the admin page run:

    $ python metab_admin.py $CONFIG_PATH

## Testing

To run a single test, run:

    $ python -m unittest tests/<desired_test>

To run all the tests, run:
    
    $ python -m unittest

If you add additional tests, the filename should begin with 'test'.
