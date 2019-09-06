Metab Import
============

## Prerequisites

Setup a Python Virtual Environment, then install the required dependencies.

    $ python3 -m venv venv
    $ source venv/bin/activate
    $ pip install -r requirements.txt


## Preparing your config file

To run the importer, you will need a config file containing two sections: one related to the Postgres database, and one related to accessing VIVO. The Postgres section includes credentials for accessing the database, the open port you used to set up the SSH tunnel, and the name of the database. The VIVO section includes credentials for accessing the VIVO query API, as well as the query endpoint and the namespace used for VIVO URIs. There is an example config file with the required fields.


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

This will produce up to four files: projects.rdf, studies.rdf, datasets.rdf, and people.rdf. These files contain the triples for each respective class. They can be uploaded to VIVO via the Add/Remove RDF data option under the Site Admin menu.

## Running the Admin Page

To start the admin page run:

    $ python metab_admin.py $CONFIG_PATH
