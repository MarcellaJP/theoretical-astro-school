FROM python:3.6-slim-buster

LABEL maintainer="Marcella Wijngaarden <marcellawijngaarden@hotmail.com>"

# The enviroment variable ensures that the python output is set straight
# to the terminal with out buffering it first
ENV PYTHONUNBUFFERED 1

WORKDIR "/ata-school"

# Install required libraries and binaries here
RUN set -ex && apt-get update && apt-get install -y --no-install-recommends \
    emacs vim nano \
    curl wget \
    git \
    make cmake \
    gcc g++ gfortran \
    # Needed for OpenMPI
    openssh-server openssh-client \
    mpich libmpich12 libmpich-dev mpich-doc \
    liblapack3 liblapack-dev libblas3 libblas-dev \
    xorg-dev perl-modules \
    libgsl-dev \
    ffmpeg && \
    echo -e "\n\nSUCCESS: done installing all required system dependencies."

# Adding the requirements ensures that Docker's layer cache kicks in :-)
ADD requirements.txt /ata-school/requirements.txt

# We also install our Python dependencies
RUN set -ex \
    pip install --upgrade pip && \
    pip install -r /ata-school/requirements.txt && \
\
    echo -e "\n\nSUCCESS: done installing all required Python dependencies."

# Create folder for software
Run set -ex && \
    mkdir /sw/headas/heasoft-6.26.1

# Install HARMPI
RUN set -ex && \
    git clone https://github.com/atchekho/harmpi /sw/harmpi && \
    cd /sw/harmpi && \
\
    # b/c to prevent future build from upstream changes
    git checkout aca57602bcff83393aeacea91c998c1adf1bcc94 && \
\
    # on this system, we need the -lm -lmpi flags - see HARMPI docs ..
    # .. so this is a search-and-replace in-place within the makefile to uncomment 'm
    sed -i -e 's/EXTRALIBS= #-lm #-lmpi/EXTRALIBS= -lm -lmpi/g' makefile && \
\
    make -j8 && \
    echo -e "\n\nSUCCESS: done installing HARMPI."

# Install Heasoft etc. for ISIS. This adds ~15 minutes to
# building the Docker image.
# If you don't intend to use ISIS, leave this uncommented.
# If you'd like to use ISIS, uncomment these lines (L60-64):
#RUN set -ex && \
#    echo -e "\n\nInstalling ISIS and prerequisites. This will take a while." && \
#    cd /install/install_isis && \
#    make -j8 software && \
#    echo -e "\n\nSUCCESS: done installing ISIS and prerequisites."

HEALTHCHECK --interval=2m --timeout=3s \
    CMD curl --fail http://localhost:8885 || exit 1

# Alternatively to CMD w/ "--allow-root", we could change the USER of container
# USER marcella

# ENTRYPOINT ["/ata-school/entrypoint.sh"]
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8885", "--allow-root"]

# Expose the port of the Docker container
EXPOSE 8885