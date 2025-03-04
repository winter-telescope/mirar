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

Installation
------------

Step 1: Install Docker

Docker is a platform for developing, shipping, and running applications in containers. You can download Docker from the following link: https://docs.docker.com/get-docker/`

Step 2: Install Docker Compose

Docker Compose is a tool for defining and running multi-container Docker applications. You can download Docker Compose from the following link: https://docs.docker.com/compose/install/

Step 3: Ensure that Docker is running

You can check if Docker is running by running the following command:

```
docker --version
```

If not, you should start Docker.

Step 4: Clone the mirar repository

You can clone the mirar repository by running the following command:

```
git clone git@github.com:winter-telescope/mirar.git
```

Step 4b (Optional): Build the Docker images

You can build a docker image, using your local copy of the mirar repository, by running the following command:

```
docker build --progress=plain --platform=linux/amd64 -t winter-telescope/mirar:latest .
```

If you want to edit the code, you should build the image afterwards. Expect building to take a few minutes. If you only want to run the latest version of the code, you can skip this step. Instead, the latest version of the pipeline will be pulled from DockerHub.

Step 5: Create a .env file

You must create a `.env` file by copying the `env.example` file in the mirar repository.

At a minimum, you will need to set the following environment variables:
* RAW_DATA_DIR: The location of the data on your machine. This is where the pipeline will read the input files.
* OUTPUT_DATA_DIR: The location of the data on your machine. This is where the pipeline will write the output files. It can be the same as RAW_DATA_DIR.

Other environment variables that are required will depend on the specific pipeline. For example, many pipelines use Astromertry.net, and will therefore require the ANET_INDEX_DIR environment variable to be set.

Step 6: Run mirar via docker-compose

You can run the docker-compose up command by running the following command:

```
docker-compose run --rm mirar-pipeline -p winter -n 20250116 -c log
```

This will remove the container after execution.

Step 6b (Optional): Run mirar via docker-compose with postgres integration

Some pipelines have integration with a postgres database. To use them, you must first specify the postgres admin credentials using your .env file.

Once this is done, you can instead use the following command to run mirar and spin up a database in the background:

```
docker-compose --profile database run --rm mirar-pipeline -p wasp -n 20250116 -c log
```

Though the main container will be removed after execution, the database will persist.

The database is configured to be accessible on your local machine via port 5433. You can
You can bring the database up and down manually with the following commands:

```
docker-compose up -d
```

and

```
docker-compose down
```
