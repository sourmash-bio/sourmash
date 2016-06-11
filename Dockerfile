FROM andrewosh/binder-base

USER root

ENV PACKAGES python-dev zlib1g git python-setuptools g++ make ca-certificates

### don't modify things below here for version updates etc.

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN pip install -U setuptools pip
RUN pip install sourmash

USER main
