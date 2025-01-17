# Use an official Python base image
ARG PSQL_MAJOR_VERSION=14
ARG PSQL_MINOR_VERSION=5
FROM postgres:${PSQL_MAJOR_VERSION}.${PSQL_MINOR_VERSION} as base

# Install required tools and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl bzip2 build-essential git ca-certificates && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
	apt-get install -y --no-install-recommends \
		# runtime requirement for using spatialite with sqlite_fdw
		libsqlite3-mod-spatialite \
		pgagent \
		postgresql-$PG_MAJOR-q3c && \
	if [ "$PG_MAJOR" -ge 14 ]; then \
		apt-get install -y --no-install-recommends postgresql-$PG_MAJOR-pgfaceting; \
	fi && \
	apt-get purge -y --auto-remove && \
	rm -rf /var/lib/apt/lists/*
