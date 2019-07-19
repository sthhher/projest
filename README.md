# Projest

Made by [Esther Vendrell]

Last updated: July 2019

## Introduction
Projest is my first program. Consist in introduce the name of a protein and obtain information about that as the most similar protein, the differents models it has or the atoms. All the information you can obtain comes from: https://www.rcsb.org/

## Description
The main file is protein.py. In spite of main functions are in extract_atoms_information.py - load protein pdb or load protein fasta - the rest of functions are in protein. 
In this way we will have functions directly linked to protein like:
Protein("2ki5").get_similar_protein() -> We will obtain fasta file, pdb file, all the information of protein in mongodb (in the init) and the similar protein of each chain. 

## Functions
There are 7 functions:
  - get_sequence_aminoacids: in Protein. Returns a list with chain and the sequence of aminoacids (1 letter). From fasta file.
  - get_chain_list: in Protein. Return a list with the differents chains there are.
  - get_aminoacid_list: in Protein. Return a list with the differents aminoacids (3 letters) with the differents sequence numbers classified in chains
  - get_similar_protein: in Protein. Return a dictionary with differents chains and the most similar protein of that chain of the protein.
  - general_dictionary: in Protein. Any return. It inserts in mongodb all the information classified. It is not callable.
  - load_protein_pdb: in extract_atoms_information. Returns the pdb file of the protein. Download the pdb file of the protein in the same directory of the files. 
  - load_protein_fasta: in extract_atoms_information. Returns the fasta file of the protein.Download the fasta file of the protein in the same directory of the files. 
