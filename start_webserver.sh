#!/usr/bin/env bash

# Extract biological data from packed archive
[ -d 'bio_data' ] && tar xfz bio_data.tar.gz

# Set default backend for Matplotlib
mkdir -p ~/.matplotlib
echo "backend: Agg" > ~/.matplotlib/matplotlibrc
# Setup ImageMagick
sed -re 's/<policy domain="coder" rights="none" pattern="PS" \/>/<!--<policy domain="coder" rights="none" pattern="PS" \/>-->/' -i /etc/ImageMagick-6/policy.xml

# Setup PostgreSQL
export POSTGRES_VERSION=12
# Accept incoming connections from all IPs ('*')
sed -re "s/#listen_addresses = 'localhost'/listen_addresses = '*'/" -i /etc/postgresql/$POSTGRES_VERSION/main/postgresql.conf
# Disable default permissions
sed -re "s/local   all             postgres                                peer/#local   all             postgres                                peer/" -i /etc/postgresql/$POSTGRES_VERSION/main/pg_hba.conf
sed -re "s/local   all             all                                     peer/#local   all             all                                     peer/" -i /etc/postgresql/$POSTGRES_VERSION/main/pg_hba.conf
sed -re "s/host    all             all             127.0.0.1\/32            md5/#host    all             all             127.0.0.1\/32            md5/" -i /etc/postgresql/$POSTGRES_VERSION/main/pg_hba.conf
sed -re "s/host    all             all             ::1\/128                 md5/#host    all             all             ::1\/128                 md5/" -i /etc/postgresql/$POSTGRES_VERSION/main/pg_hba.conf
# Permissions: Allow local connections to the database without password
echo "# Permissions: Allow local connections to the database without password" | tee -a /etc/postgresql/$POSTGRES_VERSION/main/pg_hba.conf
echo "local   all             postgres                                trust" | tee -a /etc/postgresql/$POSTGRES_VERSION/main/pg_hba.conf
echo "local   all             all                                     trust" | tee -a /etc/postgresql/$POSTGRES_VERSION/main/pg_hba.conf
echo "host    all             all             127.0.0.1/32            trust" | tee -a /etc/postgresql/$POSTGRES_VERSION/main/pg_hba.conf
echo "host    all             all             127.0.0.1/32            trust" | tee -a /etc/postgresql/$POSTGRES_VERSION/main/pg_hba.conf
# Start PostgreSQL database
service postgresql restart
# Import data
createdb -U postgres HeterologousProteinExpression
gunzip -c database/HeterologousProteinExpression.dump.sql.gz | psql -U postgres HeterologousProteinExpression

# Start Redis server
redis-server redis.conf

# Start Flask webserver
venv/bin/python web/app.py
