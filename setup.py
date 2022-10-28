import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="winterdrp",
    version="0.4.2",
    author="Viraj Karambelkar and Robert Stein",
    author_email="rdstein@caltech.edu",
    description="Data reduction pipeline for WINTER",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    keywords="astronomy image WINTER",
    url="https://github.com/winter-telescope/winterdrp",
    packages=setuptools.find_packages(),
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires='>=3.11',
    install_requires=[
        "astroplan",
        "astropy[recommended]",
        "astroquery",
        "avro-python3~=1.10.1",
        "coveralls",
        "confluent_kafka",
        "docker",
        "ephem",
        "fastavro",
        "jupyter",
        "matplotlib",
        "numpy~=1.22.2",
        "pandas",
        "penquins",
        "photutils",
        "psycopg[binary]",
        "pyfftw",
        "setuptools",
        "watchdog",
        "wget"
    ]
)