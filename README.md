Odysseus - Heterologous Gene Expression
=======================================

Odysseus is a new gene design software to optimize protein sequence data (CDS) for heterologous expression applications in common model organisms. These typical genes, adapted to the genomic characteristics of a host organism, help regulate protein expression   in such a way that a controlled decrease or increase in the rate of expression can be achieved. The software is written in Python 2.7 and realized using a probabilistic model in form of a Markov chain. The central Markov model is highly configurable by selecting pre-trained host profiles and allows to generate, for a given protein or gene sequence, sets of new synonymous DNA sequences adapted to a preset host organism. For detailed information on the functionality of the software, please see the publication linked below.

> Simm D., Popova B., Braus G. H., Waack S. and Kollmar M. (2021) Design of typical genes for heterologous gene expression.
> Scientific Reports. 12(9625). [doi: 10.1038/s41598-022-13089-1](https://doi.org/10.1038/s41598-022-13089-1)


Online access
-------------
The Odysseus gene design software can be accessed and tried out as a web-application hosted under https://odysseus.motorprotein.de.

Quick Package Installation
--------------------------

For a local installation, the user has the choice between an automated Docker installation (**recommended**) or a manual system installation.

**Docker installation**: Requires an already installed Docker environment. Run the setup script (`docker-setup.sh`) in the root directory of the cloned project. It automatically builds a Docker container and starts the Odysseus web-server running locally under http://127.0.0.1:2020: 

    $ git clone https://github.com/dsimm/Odysseus
    $ cd Odysseus
    $ ./docker-setup.sh

**Manual system installation**: The installation has been tested on Ubuntu and MacOS (HomeBrew). Below the main installation steps  are listed to install the virtual environment with the required Python modules and the Odysseus package:

    ## Setup Python project
    $ git clone https://github.com/dsimm/Odysseus
    $ cd Odysseus
    $ pip install virtualenv
    $ virtualenv venv
    $ . venv/bin/activate
    $ pip install --trusted-host pypi.python.org -r requirements.txt
    $ venv/bin/python2.7 setup.py develop

    ## System calls: Web-server and databases
    ## Start PostgreSQL database (port 5432) and Redis server (port 6379) 
    $ service postgresql restart OR pg_ctl -D /usr/local/var/postgres start -l postgres.logfile
    $ redis-server redis.conf
    
    ## Start web-server on local port 2020
    ## Run in project root directory (with virtualenv installed and activated)
    $ python web/app.py   

For detailed installation instructions on Ubuntu, see the Dockerfile and the further software requirements listed below.

System Requirements
-------------------

As a working environment for the application the following software packages need to be installed on the system:

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
