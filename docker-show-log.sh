#!/usr/bin/env bash

echo "" && echo "Container: odysseus"
docker exec -it odysseus sh -c "tail -n 100 -f log/*.log"
