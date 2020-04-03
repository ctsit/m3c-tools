Metab Import
============

## Prerequisites

Setup a Python Virtual Environment, then install the required dependencies.

    $ python3 -m venv venv
    $ source venv/bin/activate
    $ pip install -r requirements.txt

Also obtain a PubMed API token for increased API limits. See [this site](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) for more information.

## Preparing your config file

To run the importer, you will need a config file containing two sections: one
related to the Postgres database, and one related to accessing VIVO. The
Postgres section includes credentials for accessing the database, the open port
you used to set up the SSH tunnel, and the name of the database. The VIVO
section includes credentials for accessing the VIVO query API, as well as the
query endpoint and the namespace used for VIVO URIs. There is an example config
file with the required fields. Download a CSV export from the
[Metabolomics Software Tools](https://docs.google.com/spreadsheets/d/1a096jlzbAwTxUdvTtJB2bwfkBRdDcMKrkywS8YZf344/edit?usp=sharing)
and note the file name in the config for `csv_tools`.


## Connect to the database

The first step to using the importer is setting up an SSH tunnel to
stage.vivo.metabolomics.info to access the Postgres database with the
information from Metabolomics Workbench. You will need to set up the tunnel in
its own terminal tab.

    $ ssh -N -L 5432:localhost:5432 stage.vivo.metabolomics.info


## Running the Pre-fill script

Once the tunnel is prepared, you can run the pre-fill script.

    $ python m3c/prefill.py $CONFIG_PATH

This pre-fills the supplemental tables with necessary information like people
and organizations.


## Running the Importer

Next, run:

    $ python metab_import.py $CONFIG_PATH

This will produce up to four files: projects.rdf, studies.rdf, datasets.rdf,
and people.rdf. These files contain the triples for each respective class.


## Running the Publication Fetcher

The Publication Fetcher tries to find all authors' publications by using
Harvard's Catalyst service or PubMed. Then downloads the XML summaries from
PubMed and adds them to the supplemental database. Use the admin page to
mark publications for inclusion and exclusion. (At least one PMID and an
affiliation is required to use Catalyst).

    $ python pubfetch.py $CONFIG_PATH


## Running the Admin Page

To start the admin page run:

    $ python metab_admin.py $CONFIG_PATH


## Testing

To run all the tests, run:

    $ python -m unittest

To run a single test, run:

    $ python -m unittest tests/<desired_test>

If you add additional tests, the filename should begin with 'test'.

