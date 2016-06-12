# for binder
FROM andrewosh/binder-base

USER root

ENV PACKAGES python-dev zlib1g git python-setuptools g++ make ca-certificates
ENV SOURMASH_VERSION master

### don't modify things below here for version updates etc.

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN cd /home && \
    git clone https://github.com/dib-lab/sourmash.git -b ${SOURMASH_VERSION} && \
    cd sourmash && \
    python setup.py install

USER main
