#!/usr/bin/env bash

# Extract biological data from packed archive
[ ! -d 'bio_data' ] && tar xfz bio_data.tar.gz

# Build docker image of project with composer-preferences
docker-compose build
# docker-compose -f docker-compose.yml build --no-cache

# Start image as container via composer
# docker-compose up
# docker-compose up -d
mkdir -p log
docker-compose up > log/development.log &

# Log into running docker container
# docker exec -it odysseus "/bin/bash"
