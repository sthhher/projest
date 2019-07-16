import requests
import urllib.request
from urllib.request import urlopen
from xml.dom import minidom
from Bio.SeqUtils import seq1
from Bio import SeqIO
from itertools import groupby
import string

import pymongo
from pymongo import MongoClient
client = MongoClient('localhost', 27017)
db = client["proyecto"]
collection = db["proteinas"]

import extract_atoms_information
from chain import Chain
from aminoacid import Aminoacid
from atom import Atom
from model import Model

#SETDEFAULT

class Protein(object):

    def __init__(self, protein_name):
        self.protein_name = protein_name
        extract_atoms_information.load_protein_pdb(self.protein_name) #It creates automatically so, we dont have to do it.
        extract_atoms_information.load_protein_fasta(self.protein_name)
    
    def get_sequence_aminoacids (self): #NO ACABADA
        record = []
        record = list(SeqIO.read(self.protein_name + "_fasta", "fasta"))
        print(record[0].seq[0:])

    def chain_list(self): #Return the differents chains there are
        list_differents_chains = []
        with open(self.protein_name, 'r+') as text_file:
            for line in text_file:
                if line.startswith("ATOM"):
                    Chain.chain_identifier = line[21]
                    Chain.chain_identifier = Chain.chain_identifier.strip()
                    if Chain.chain_identifier not in list_differents_chains:
                        list_differents_chains.append(Chain.chain_identifier)
            return list_differents_chains

    def aminoacid_list(self):
        with open(self.protein_name, 'r+') as text_file:
            for line in text_file:
                if line.startswith("ATOM"):
                    Aminoacid.residue_name = line[17:20]
                    residue_name = Aminoacid.residue_name.strip()

                    Aminoacid.residue_sequence_number = line[22:26]
                    residue_sequence_number = Aminoacid.residue_sequence_number.strip()
                    return residue_name, residue_sequence_number

    def get_similar_protein(self): #VA BIEN -> HACER UN TEST Y UN IF EN LA FUNCION PARA SI UNA CADENA DE PROTEINA ESTA MAL
        dicc_hit_def = {}
        open_fasta = list(SeqIO.parse(self.protein_name + "_fasta", 'fasta'))
        chain = (len(open_fasta)-1)
        while chain >=0:
            fasta_sequence = open_fasta[chain].seq
            url = "https://www.rcsb.org/pdb/rest/getBlastPDB1?sequence=" + str(fasta_sequence) + "&eCutOff=10.0&matrix=BLOSUM62&outputFormat=XML" #TENGO QUE HACER QUE ME DIGA EL NOMBRE DE LA PROTEINA
            request = requests.get(url) #Here is where im getting the status_code
            file_protein_pdb = "similar_protein"

            if request.status_code == 200: #Check out if a file exists
                Download = urlopen(url) #Download pdb
                file = open(file_protein_pdb, "wb") #I want to write in bytes in file_protein_pdb
                file.write(Download.read()) 
                file.close()
                doc = minidom.parse(file_protein_pdb) #This is for read the xml
                hits = doc.getElementsByTagName("Hit")
                hit_def = hits[0].getElementsByTagName("Hit_def")[0] #Define Hit_def
                hit_def = hit_def.firstChild.data #I want values of hit_def
                hit_def = hit_def[0:4] #I only want 4 first values of hit_def
                dicc_hit_def[open_fasta[chain].id[5]]=hit_def #Add to the dictionary
                chain = chain-1
        return dicc_hit_def

    def create_model_dictionary(self):
        model_count = 0
        model_dictionary = {}
        with open(self.protein_name, 'r+') as text_file:
            for line in text_file:
                if line.startswith("MODEL"):
                    #MODEL
                    model_count +=1

                    Model.model_identifier = line[11:17]
                    model_identifier = Model.model_identifier.strip()

                    model_dictionary.setdefault(model_identifier,{})
        if model_count != 0:
            return model_dictionary
        else:
            model_dictionary = {"1":{}}
            return model_dictionary

    def general_dictionary(self):
        value_dictionary = 0

        protein = {"_id":self.protein_name}

        with open(self.protein_name, 'r+') as text_file:
            for line in text_file:
                if line.startswith("ATOM"):

                    value_dictionary +=1

                    #CHAIN
                    Chain.chain_identifier = line[21]
                    chain_identifier = Chain.chain_identifier.strip()

                    #AMINOACID
                    Aminoacid.residue_name = line[17:20]
                    residue_name = Aminoacid.residue_name.strip()

                    Aminoacid.residue_sequence_number = line[22:26]
                    residue_sequence_number = Aminoacid.residue_sequence_number.strip()

                    #ATOM
                    Atom.atom_name = line[12:16]
                    atom_name = Atom.atom_name.strip()

                    Atom.x_coordinate = line[30:38]
                    x_coordinate = float(Atom.x_coordinate.strip())

                    Atom.y_coordinate = line[38:46]
                    y_coordinate = float(Atom.y_coordinate.strip())

                    Atom.z_coordinate = line[46:54]
                    z_coordinate = float(Atom.z_coordinate.strip())

                    Atom.occupancy = line[54:60]
                    occupancy = float(Atom.occupancy.strip())

                    Atom.temperature_factor = line[60:66]
                    temperature_factor = float(Atom.temperature_factor.strip())

                    Atom.element_symbol = line[76:78]
                    element_symbol = Atom.element_symbol.strip()
                    
                    #SETDEFAULTs
                    chain_letter = protein.setdefault(chain_identifier, {})
                    residue_letter = chain_letter.setdefault(residue_name, {})
                    residue_number = residue_letter.setdefault(residue_sequence_number, {})
                    atom_letter = residue_number.setdefault(element_symbol, {})

                    atom_letter[atom_name] = {
                        "x_coordinate":x_coordinate, 
                        "y_coordinate":y_coordinate, 
                        "z_coordinate":z_coordinate, 
                        "occupancy":occupancy, 
                        "temperature_factor":temperature_factor,
                    }

                    residue_number[element_symbol] = atom_letter
        return protein

dictionary = Protein("2ki5").general_dictionary()
collection.insert_one(dictionary)
print(dictionary)
