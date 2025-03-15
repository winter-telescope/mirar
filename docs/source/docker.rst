Docker Integration
============


Installing the package
----------------------

The simplest way to run mirar is via Docker.
Images for various versions are available on Docker Hub.
You can run the mirar code within a Docker container by following the steps below.
You first specify the directories on your local machine where the data is stored and where the output will be written, using a `.env` file.
You can then choose the pipeline to run, the night to process, and the configuration to use as command line arguments.


Step 1: Install Docker
......................

Docker is a platform for developing, shipping, and running applications in containers. You can download Docker from the following link: https://docs.docker.com/get-docker/`

Step 2: Install Docker Compose
...............................

Docker Compose is a tool for defining and running multi-container Docker applications. You can download Docker Compose from the following link: https://docs.docker.com/compose/install/

Step 3: Ensure that Docker is running
.....................................

You can check if Docker is running by running the following command:

.. code-block:: bash

    docker --version

If not, you should start Docker.

Step 4: Clone the mirar repository
..................................

You can clone the mirar repository by running the following command:

.. code-block:: bash

    git clone git@github.com:winter-telescope/mirar.git

Step 4b (Optional): Build the Docker image
............................................

If you like, you can build a docker image using your local copy of the mirar repository. You can do so by running the following command:

.. code-block:: bash

    docker build --progress=plain --platform=linux/amd64 -t winter-telescope/mirar:latest .

If you edit the code locally, you should build the image afterwards. Expect building to take a few minutes.

If you only want to run the latest version of the code, you can skip this step. Instead, the latest version of the pipeline will be pulled from DockerHub.

Step 5: Create a .env file
...........................

You must create a `.env` file by copying the `env.example` file in the mirar repository.

At a minimum, you will need to set the following environment variables:
- RAW_DATA_DIR: The location of raw data on your machine. This is where the pipeline will read the input files.
- OUTPUT_DATA_DIR: The location for mirar to output data on your machine. It can be the same as RAW_DATA_DIR.

Other environment variables that are required will depend on the specific pipeline.
For example, many pipelines use Astromertry.net, and will therefore require the ANET_INDEX_DIR environment variable to be set.
You will need to follow the instructions in :doc:`installation` for downloading these files.

Step 6: Run mirar via docker-compose
.....................................

You can run the docker-compose up command by running the following command:

.. code-block:: bash

    docker-compose run --rm mirar-pipeline -p wasp -n 20250116 -c log

This will reduce data under the subdirectory /wasp/20250116/raw of RAW_DATA_DIR, using the 'log' configuration.
Files will be output in a subdirectory under /wasp/20250116 of OUTPUT_DATA_DIR.
Once processing is completed, the container will close.

Step 6b (Optional): Run mirar via docker-compose with postgres integration
...........................................................................

Some pipelines have integration with a postgres database, for example WINTER.
To use them, you must first be sure to specify the postgres admin credentials using your .env file.

Once this is done, you can use the a slightly different command to run mirar and spin up a database in the background:

.. code-block:: bash

    docker-compose --profile database run --rm mirar-pipeline -p winter -n 20250116 -c log

This command will in general do everything described in step 6, but will also spin up a postgres database in the background.
Though the main container will be removed after execution, the database will persist.

The database is configured to be accessible on your local machine via port 5433. You can interact with it using a tool like pgAdmin or DBeaver, or a local psql client.

You can bring the database up and down manually with the following commands:

.. code-block:: bash

    docker-compose --profile database up -d

and

.. code-block:: bash

    docker-compose --profile database down

Bringing the database up and down will not affect the data in the database.
