################################################################################
# Install basic system for Python2.7
FROM	ubuntu:20.04

# Optionally set a maintainer name to let people know who made this image.
MAINTAINER Dominic Simm <dominic.simm@cs.uni-goettingen.de>

# Install dependencies:
# - build-essential: Compiler environment gcc
# - python-dev: Compilation dep. Python module biopython
# - pkg-config: Compilation dep. Python module matplotlib
# - libfreetype6-dev: Compilation dep.: Python module matplotlib  (dev. files)
# - libpng-dev: Compilation dep.: Python module matplotlib (dev. files)
# - libpq-dev: Compilation dep.: Python module psycopg2 (dev. files)
# - gnuplot: GNUplot
# - graphviz: Graph lib
# - imagemagick: Image library (convert)
ARG		DEBIAN_FRONTEND=noninteractive
RUN		apt-get update --fix-missing && \
		apt-get install -y apt-utils build-essential coreutils pkg-config
RUN		apt-get install -y python-dev libfreetype6-dev libpng-dev libpq-dev
RUN		apt-get install -y curl wget rsync git sudo nano screen
RUN		apt-get install -y gnuplot graphviz imagemagick
RUN		apt-get install -y postgresql redis
# Note: The official Debian and Ubuntu images automatically 'apt-get clean'
# after each 'apt-get'

# Set environment variables to store where the app is installed inside of
# the Docker image.
ENV		APP_PATH /heprex
RUN		mkdir -p $APP_PATH
ENV		VENV $APP_PATH/venv/bin
ENV		TMP_PATH /tmp

# This sets the context of where commands will be ran in and is documented
# on Docker's website extensively.
WORKDIR	$APP_PATH

# Ensure modules are cached and only get updated when they change. This will
# drastically increase build times when your modules do not change.
ADD		requirements.txt requirements.txt

# Install virtualenv with all specified modules in package_list.txt
RUN		wget https://bootstrap.pypa.io/pip/2.7/get-pip.py && python get-pip.py
RUN		python2.7 -mpip install -U pip
RUN		python2.7 -mpip install -U virtualenv
RUN		virtualenv --python=/usr/bin/python2.7 venv
RUN		$VENV/pip install -U numpy==1.10.1
RUN		$VENV/pip install --trusted-host pypi.python.org -r requirements.txt

# Copy in the application code from your work station at the current directory
# over to the working directory.
COPY	. .

# Install HePrEx module
RUN		$VENV/python setup.py develop

# Install additional software
RUN		apt-get install -y libgsl23
RUN		wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_20_04/viennarna_2.4.18-1_amd64.deb && \
        dpkg -i viennarna_2.4.18-1_amd64.deb && rm viennarna_2.4.18-1_amd64.deb

# Make port 2020 available to the world outside this container
EXPOSE	2020

# Expose a volume so that the container will be able to read files from other Docker-containers.
VOLUME	["$TMP_PATH"]

# ENTRYPOINT [ "venv/bin/python", "web/app.py"]
# CMD		nohup venv/bin/python web/app.py > /dev/null 2>&1 &
