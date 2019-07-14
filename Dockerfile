FROM python:3.6-slim

LABEL maintainer="Marcella Wijngaarden <marcellawijngaarden@hotmail.com>"

# The enviroment variable ensures that the python output is set straight
# to the terminal with out buffering it first
ENV PYTHONUNBUFFERED 1

WORKDIR "/ata-school"

# Adding the requirements ensures that Docker's layer cache kicks in :-)
ADD requirements.txt /ata-school/requirements.txt

# The container runs Debian. We install required libraries here
RUN set -ex \
    \
    && apt-get update && apt-get install -y --no-install-recommends \
    curl \
    gcc \
    gfortran \
    ffmpeg

# We also install our Python dependencies
RUN set -ex \
    pip install --upgrade pip && \
    pip install -r /ata-school/requirements.txt && \
\
    echo -e "\n\nSUCCESS: done packaging all required dependencies."

HEALTHCHECK --interval=2m --timeout=3s \
    CMD curl --fail http://localhost:1337 || exit 1

# Alternatily to CMD w/ "--allow-root", we could change the USER of container
# USER marcella

# ENTRYPOINT ["/ata-school/entrypoint.sh"]
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=1337", "--allow-root"] 

# Expose the port of the Docker container
EXPOSE 1337
