# Projest

Made by [Esther Vendrell]

Last updated: July 2019

## Introduction
Projest is my first program. Consist in introduce the name of a protein and obtain information about that as the most similar protein, the differents models it has or the atoms. All the information you can obtain comes from: https://www.rcsb.org/

## Description
The main file is protein.py. In spite of main functions are in extract_atoms_information.py - load protein pdb or load protein fasta - the rest of functions are in protein. 
In this way we will have functions directly linked to protein like:
Protein("2ki5").get_similar_protein() -> We will obtain fasta file, pdb file, all the information of protein in mongodb (in the init) and the similar protein of each chain. 


