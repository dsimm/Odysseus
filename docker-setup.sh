#!/usr/bin/env bash

# Build docker image of project with composer-preferences
docker-compose build
# docker-compose build --no-cache --force-rm

# Start image as container via composer
# docker-compose up
# docker-compose up -d
mkdir -p log
docker-compose up > log/development.log &

# Log into running docker container
# docker exec -it odysseus "/bin/bash"
