FROM python:3.6-slim

LABEL maintainer="Marcella Wijngaarden <marcellawijngaarden@hotmail.com>"

# The enviroment variable ensures that the python output is set straight
# to the terminal with out buffering it first
ENV PYTHONUNBUFFERED 1

WORKDIR "/ata-school"

# The container runs Debian. We install required libraries here
RUN set -ex \
    \
    && apt-get update && apt-get install -y --no-install-recommends \
    vim \
    curl \
    make \
    gcc \
    gfortran \
    openmpi-bin openmpi-doc libopenmpi-dev \
    git \
    ffmpeg

# Adding the requirements ensures that Docker's layer cache kicks in :-)
ADD requirements.txt /ata-school/requirements.txt

# We also install our Python dependencies
RUN set -ex \
    pip install --upgrade pip && \
    pip install -r /ata-school/requirements.txt && \
\
    echo -e "\n\nSUCCESS: done packaging all required dependencies."

# Install HARMPI
RUN set -ex && \
    git clone https://github.com/atchekho/harmpi /harmpi && \
    cd /harmpi && \
    sed -i -e 's/EXTRALIBS= #-lm #-lmpi/EXTRALIBS= -lm -lmpi/g' makefile && \
    make -j8

HEALTHCHECK --interval=2m --timeout=3s \
    CMD curl --fail http://localhost:1337 || exit 1

# Alternatily to CMD w/ "--allow-root", we could change the USER of container
# USER marcella

# ENTRYPOINT ["/ata-school/entrypoint.sh"]
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=1337", "--allow-root"] 

# Expose the port of the Docker container
EXPOSE 1337