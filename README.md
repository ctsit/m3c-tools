M3C Tools
=========

## Installation

    $ pip install git+https://github.com/ctsit/m3c-tools.git

## PubMed API Token

[Obtain a PubMed API token](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) for increased API limits.

## Configuration

To run the tools, you will need a YAML configuration file containing several
properties. See `config-example.yaml` for an example.

## Database Tunnel

If Postgres is behind a firewall, you can setup an SSH tunnel to connect to it:

    $ ssh -L 5432:localhost:5432 stage.vivo.metabolomics.info

Where `stage.vivo.metabolomics.info` is the host running Postgres and `5432` is
the standard Postgres port.


## Run the Pre-fill script

The `prefill` command pre-fills the Supplemental tables with necessary
information like people and organizations.

Run it:

    $ m3c prefill $CONFIG_PATH


## Running the Importer

Next, run:

    $ m3c generate $CONFIG_PATH

This will produce up to several N-Triples files including: projects.nt,
studies.nt, datasets.nt, and people.nt. These files contain the triples for
each respective type.


## Running the Publication Fetcher

The Publication Fetcher tries to find all authors' publications by using
Harvard's Catalyst service or PubMed. Then downloads the XML summaries from
PubMed and adds them to the supplemental database. Use the admin page to
mark publications for inclusion and exclusion. (At least one PMID and an
affiliation is required to use Catalyst).

    $ m3c pubfetch $CONFIG_PATH


## Starting the Admin Forms server

To start the Admin Froms server, run:

    $ m3c serve $CONFIG_PATH

They should be accessible at: http://localhost:5000/


## Development

Download the code, setup a virtual environment, and configure it for
development.

    $ git clone https://github.com/ctsit/m3c-tools.git
    $ cd m3c-tools/
    $ python3 -m venv venv
    $ source venv/bin/activate
    $ python setup.py develop


## Testing

To run all the tests, run:

    $ python -m unittest

To run a single test, run:

    $ python -m unittest tests/<desired_test>

If you add additional tests, the filename should begin with 'test'.

