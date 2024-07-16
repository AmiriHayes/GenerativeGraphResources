# Network Vaccination

This repository contains the code and data used in two recent publications that used mobile device data to model the spread of COVID-19 throughout a partially vaccinated population:
- [Protecting the Most Vulnerable by Vaccinating the Most Active](https://www.rand.org/pubs/perspectives/PEA1068-1.html) (RAND Perspective report)
- Modeling the Impact of Social Distancing and Targeted Vaccination on the Spread of COVID-19 through a Real City-Scale Contact Network  (forthcoming).

## Project description
Motivated by the (then) imminent availability of COVID-19 vaccines, this repository was used to compare the efficacies of different network-based vaccination strategies. It is a well-known result in network science that targeted strategies, where more connected nodes are vaccinated first, can sometimes far outperform the strategy where nodes are vaccinated uniformly at random. <sup id="pastor-satorras">[1](#f1)</sup> It also seems plausible that vaccinating nodes with low connectivities is an even worse strategy.

We investigated these hypotheses using a real-world person-to-person contact network derived from mobile device data obtained from the company Uber Media. Rather than attempting to faithfully simulate the spread of COVID-19 on this network, we instead simulated two rather simple contagion models (SIR and SEIR). Consequently, we are simulating a very idealized contagion model on a realistic network - in future work we would like to improve this by considering more realistic contagion models.

## How to use this repository

### Installation
- Download the repo: `git clone https://github.com/RANDCorporation/network_vaccination`
- Navigate to repo: `cd network_vaccination`
- Install required packages: `pip install -r requirements.txt`

It is recommended to use a conda environment, in which case the command sequence should be modified to
- `conda create -n network_env`
- `conda activate network_env`
- `git clone https://github.com/RANDCorporation/network_vaccination`
- `cd network_vaccination`
- `pip install -r requirements.txt`

### Data
The mobile device data was used to create 3 different anonymized contact networks:
- G_pre capturing interactions before social distancing measures were enacted
- G_post capturing interactions during social distancing
- G_superposition representing a combined network which is described in the [Perspective report](https://www.rand.org/pubs/perspectives/PEA1068-1.html) in more detail

The weighted adjacency lists for these anonymized contact networks are located in the `data` directory.

### Code
The repository is primarily organized into multiple Jupyter notebooks, each accomplishing a different task.
  - `Contact Network Analysis.ipynb`: This notebook analyzes the person-to-person contact networks used in the simulations.
  - Code used in [Protecting the Most Vulnerable by Vaccinating the Most Active](https://www.rand.org/pubs/perspectives/PEA1068-1.html):
    - `SIR Model.ipynb`: This notebook carries out the SIR simulation over G_superposition.
    - `SIR Data Analysis.ipynb`: This notebook analyzes the results of the simulation and makes some plots.
  - Code used in Using Mobile Device Data to Measure The Impact of Social Distancing Measures and Optimal Vaccination Strategies for COVID-19:
    - `SEIR Model.ipynb`: This notebook carries out the SEIR simulation over the G_pre and G_post contact networks.
    - `SEIR Data Analysis.ipynb`: This notebook analyzes the results of the simulation and makes some plots.

The `utils.py` folder also contains some useful functions used by the various notebooks.

### Contact
Code developed and maintained by Gavin Hartnett (hartnett@rand.org).

<b id="f1">1</b> Pastor-Satorras, Romualdo, and Alessandro Vespignani. "Immunization of complex networks."  [Physical review E 65.3 (2002): 036104](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.036104). [ArXiv link](https://arxiv.org/abs/cond-mat/0107066). [↩](#a1)
