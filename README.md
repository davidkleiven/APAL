![Build status](https://travis-ci.org/davidkleiven/APAL.svg?branch=master)

# APAL
Alloy Phasefield Abstraction Layer. The aim of APAL is to provide a simple Python interface to phase field models. 
It is solely based on the [MMSP project](https://github.com/mesoscale/mmsp). Supported equations

    * Cahn-Hilliard
    * Cahn-Hilliard-Ginzburg-Landau (aimed at precipitate modelling)

# Installation
First, install the dependencies (scipy and cython)

```bash
pip install -r requirements.txt --user
```
then compile and install the APAL package

```bash
pip install . --user
```
