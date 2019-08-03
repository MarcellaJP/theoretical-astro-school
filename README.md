# My repo for the Advancing Theoretical Astrophysics School

Repository for code related to the [Advancing Theoretical Astrophysics summer school](https://collectiveastronomy.github.io/advancingtheoastro/) (15 - 26th July 2019) at the University of Amsterdam.

## Setup

I am using a Docker container to include the full environment (from OS to Python dependencies to code) used for this school. My host computer OS is Windows 10 and in the container I am running Debian Linux. I chose to use a Docker setup because I am experimenting with using and publishing Docker containers for completely reproducible science environments.

This repository includes the tutorials used in the Advancing Theoretical Astrophysics summer school and the simulation output data necessary to run the tutorial notebooks. You can also rerun simulations in the Docker container to create new data to use in the tutorials (or make movies!). You can run Harmpi and ISIS in the container (see install notes).

## Install & use instructions

To run this Docker environment (and all the code in it), you only need to install the latest versions of `Docker` and `docker-compose`. 

To start the Docker container (incl. serving Jupyter notebook at  `localhost:8885`), run the following command in the top level of this repository:

`docker-compose -f docker-compose.yml up --build`

And to start the container shell:

`docker exec -it ata-school-container bash`

See [the Docker documentation](https://docs.docker.com/) and [commonly used commands](https://towardsdatascience.com/15-docker-commands-you-should-know-970ea5203421).


## Quick links

* [Schedule](https://github.com/collectiveastronomy/advancingtheoastro/blob/master/ATABlockSchedule.pdf)
* [Course curriculum](https://github.com/collectiveastronomy/advancingtheoastro/blob/master/ata_curriculum.pdf)
* [Course materials](https://github.com/collectiveastronomy/ATAMaterials)

