Docker Integration
============


Installing the package
----------------------

The simplest way to run mirar is via Docker. Images for various versions are available on Docker Hub.

The basic set up is the following:

* A docker image is created with the mirar package installed.
* A second docker image is created with a postgres database installed, alongside q3c.
* As a user, you can copy the docker-compose.yaml file from the mirar repository. This file ensures that two images play nicely together.
* Important environment variables are set via a .env file. These include the database name, user, password, and the location of the data on your machine.
* After configuring a .env file, you can run the docker-compose up command. This will start the two images and create a network between them.

Important environment variables that must be set:
* PG_ADMIN_PWD: The password for the postgres admin user. Even if you only use a pipeline without psql inegration, you must set this to some value.
* RAW_DATA_DIR: The location of the data on your machine. This is where the pipeline will read the input files.
* OUTPUT_DATA_DIR: The location of the data on your machine. This is where the pipeline will write the output files.

The following environment variables are often needed:
* ANET_INDEX_DIR: The location of the anet index on your machine. This is where the pipeline will read the anet index files.


Commands to run
----------------------

To build or rebuild the images, run the following command:

```
docker build --progress=plain --platform=linux/amd64 -t winter-telescope/mirar:v0.20.0 .
```

and:

```
docker build --progress=plain --platform=linux/amd64 -t winter-telescope/mirar-db:v0.20.0 .
```

To reduce data:

```
docker-compose run --rm mirar-pipeline -p winter -n 20250116 -c log
```

You do not need to explicitly start up the database, since the `docker-compose run` command will start it up for you.
This will remove the container after execution. However, the database will persist.

The database is configured to be accessible on your machine via port 5433.
You can bring the database up and down manually with the following commands:

```
docker-compose up -d
```

and

```
docker-compose down
```
