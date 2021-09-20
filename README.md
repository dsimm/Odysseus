Odysseus - Heterologous Protein Expression
==========================================

Odysseus is a new gene design software to adapt protein sequence data (CDS) for heterologous expression applications in common model organisms. The software is written in Python 2.7 and realized using a probabilistic model in form of a Markov chain as DNA sequence generator producing adapted synonymous gene sequences of the input protein or gene. 

Usage
-----

System call: web-server on local port 2020

    # Run in project root directory (with virtualenv installed and actived)
    $ python web/app.py


Quick Package Installation
--------------------------

The user has the choice between a Docker installation (recommended) or a manual system installation.

**Docker installation**: Run the following script (docker-setup.sh) in the root directory of the project. It 
automatically builds a docker image and starts the Odysseus local web-server running under http://127.0.0.1:2020: 

    $ git clone https://github.com/dsimm/Odysseus
    $ cd odysseus
    $ ./docker-setup.sh

Manual system installation (for detailed instructions see the Dockerfile):

    # Setup Python project
    $ git clone https://github.com/dsimm/Odysseus
    $ cd odysseus
    $ pip install virtualenv
    $ virtualenv venv
    $ . venv/bin/activate
    $ pip install --trusted-host pypi.python.org -r requirements.txt
    $ venv/bin/python2.7 setup.py develop

    # Start PostgreSQL database (port 5432) and Redis Server (port 6379) 
    $ pg_ctl -D /usr/local/var/postgres start -l postgres.logfile
    $ redis-server redis.conf


System Requirements
-------------------

As a working environment for the web-application the following software packages need to be installed on the system:

* Python2.7 with pip
* PostgreSQL (>= 11.4)
* Redis (>= 5.0.x)
* gnuplot (>= 5.2)
* graphviz (>= 2.4)
* imagemagick (>= 7.0.x)
* ViennaRNA (>= 2.4.4.1)


Detailed Installation Instructions
----------------------------------

Coming soon ...
