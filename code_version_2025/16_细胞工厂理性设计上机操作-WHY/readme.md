# Hands-on Practice in Rational Strain Design
## 1. Introduction
First, demonstrate how to conduct rational strain design using cellular digital models such as GEMs and enzyme-constrained metabolic models with the Cobrapy toolkit and Python programming language. Under the guidance of instructors and teaching assistants, each student selects a specific bio-based chemical compound and a microbial chassis. They then use GEMs or enzyme-constrained metabolic models in practical computer operations to rationally design the microbial chassis strains, predicting gene modification targets required to enhance the production of the target compound.
## 2. Setup Environment
### 1. Create a conda environment
```bash
conda create -n strain_design
conda activate strain_design
```
### 2. Install required packages
```bash
conda install -c conda-forge cobra, seaborn, pandas
```
## 3.Tutorials for rantional strain design
### 1. [Strain design with GEMs by FSEOF algorithm](./strainDesign_with_ecGEMs.py)
### 2. [Strain design with GEMs by optForce algorithm](./optForce_with_GEM.py)
### 3. [Strain design with enzyme-constrained GEMs by ecFSEOF algorithm](./strainDesign_with_ecGEMs.py)



