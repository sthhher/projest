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

class Protein(object):

    def __init__(self, protein_name):
        self.protein_name = protein_name
        extract_atoms_information.load_protein_pdb(self.protein_name) #Que se cree automaticamente para no tener que hacerlo en cada funcion
        extract_atoms_information.load_protein_fasta(self.protein_name)
        #extract_atoms_information.insert_proteins_mongo(self.protein_name)
    
    def get_sequence_aminoacids (self): #NO ACABADA
        record = []
        record = list(SeqIO.read(self.protein_name + "_fasta", "fasta"))
        print(record[0].seq[0:])

    def chain_list(self): #Return the differents chains there are
        list_differents_chains = []
        cadenas = {"_id":self.protein_name}
        with open(self.protein_name, 'r+') as text_file:
            for line in text_file:
                if line.startswith("ATOM"):
                    Chain.chain_identifier = line[21]
                    Chain.chain_identifier = Chain.chain_identifier.strip()
                    if Chain.chain_identifier not in list_differents_chains:
                        list_differents_chains.append(Chain.chain_identifier)
                        cadenas[Chain.chain_identifier] = {}
            return cadenas

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

    def general_dictionary(self):
        value_dictionary = 0
        residue_A = {}
        residue_B = {}
        residue_C = {}
        residue_D= {}
        residue_E = {}
        residue_F = {}
        residue_G = {}
        residue_H = {}
        residue_I = {}
        residue_J = {}
        residue_K = {}
        residue_L = {}
        residue_M = {}
        residue_N = {}
        residue_O = {}
        residue_P= {}
        residue_Q = {}
        residue_R = {}
        residue_S = {}
        residue_T = {}
        residue_U = {}
        residue_V = {}
        residue_W = {}
        residue_X = {}
        residue_Y = {}
        residue_Z = {}

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
                    Atom.atom_serial_number = line[6:11]
                    atom_serial_number = Atom.atom_serial_number.strip()

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
                    
                    if chain_identifier == "A":
                        residue_A[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_A

                    elif chain_identifier == "B":
                        residue_B[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_B

                    elif chain_identifier == "C":
                        residue_C[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_C

                    if chain_identifier == "D":
                        residue_D[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_D
    
                    elif chain_identifier == "E":
                        residue_E[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_E

                    elif chain_identifier == "F":
                        residue_F[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                    elif chain_identifier == "G":
                        residue_G[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_G

                    elif chain_identifier == "H":
                        residue_H[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_H

                    elif chain_identifier == "I":
                        residue_I[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_I
    
                    elif chain_identifier == "J":
                        residue_J[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_J

                    elif chain_identifier == "K":
                        residue_K[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_K

                    elif chain_identifier == "L":
                        residue_L[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_L

                    elif chain_identifier == "M":
                        residue_M[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_M

                    elif chain_identifier == "N":
                        residue_N[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_N
    
                    elif chain_identifier == "O":
                        residue_O[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_O

                    elif chain_identifier == "P":
                        residue_P[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_P

                    elif chain_identifier == "Q":
                        residue_Q[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_Q

                    elif chain_identifier == "R":
                        residue_R[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_R

                    elif chain_identifier == "S":
                        residue_S[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_S
    
                    elif chain_identifier == "T":
                        residue_T[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_T

                    elif chain_identifier == "U":
                        residue_U[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_U

                    elif chain_identifier == "V":
                        residue_V[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_V

                    elif chain_identifier == "W":
                        residue_W[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_W

                    elif chain_identifier == "X":
                        residue_X[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_X

                    elif chain_identifier == "Y":
                        residue_Y[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_Y

                    elif chain_identifier == "Z":
                        residue_Z[residue_name + " - " + str(value_dictionary)] = {
                            "chain_identifier":chain_identifier,
                            "residue_sequence_number":residue_sequence_number,
                            "residue_name":residue_name,
                            "ATOM":{
                                "atom_serial_number":atom_serial_number, 
                                "atom_name":atom_name, 
                                "x_coordinate":x_coordinate, 
                                "y_coordinate":y_coordinate, 
                                "z_coordinate":z_coordinate, 
                                "occupancy":occupancy, 
                                "temperature_factor":temperature_factor,
                                "element_symbol":element_symbol
                            }
                        }

                        protein[chain_identifier] = residue_Z
        return protein


dictionary = Protein("4M8B").general_dictionary()
collection.insert_one(dictionary)
print(dictionary)
