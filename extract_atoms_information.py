import urllib.request
from urllib.request import urlopen
import requests
import pymongo


def load_protein_pdb(name_protein):
    fileDownload = "https://files.rcsb.org/download/" + name_protein + ".pdb" #Get the link
    request = requests.get(fileDownload) #Here is where im getting the status_code
    
    if request.status_code == 200: #Check out if a file exists
        Download = urlopen(fileDownload) #Open url
        file = open(name_protein, "wb") #Open the file, write binary
        file.write(Download.read())

def load_protein_fasta(name_protein):
    url = "https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=" + name_protein + "&compressionType=uncompressed"
    request = requests.get(url)

    if request.status_code == 200: #Check out if a file exists
        Download = urlopen(url) #Open url
        file = open(name_protein + '_fasta', "wb") #Open the file, write binary
        file.write(Download.read()) 
    

def insert_proteins_mongo(protein_name):
        myclient = pymongo.MongoClient("mongodb://localhost:27017/")
        mydb = myclient["proyecto"]
        mycol = mydb["proteinas"]
        mylist = [
            {"_id":protein_name,
            protein_name:
                {"Chain_A":
                    {"1_MET":
                        {"Atom_1":
                            {"atom_name":"AA", 
                             "x_coordinate":3.455
                            },
                        "Atom_2":
                            {"atom_name":"T"
                            }
                        },
                    "2_TYR":
                        {"Atom_1":
                            {"atom_name":"T"
                            }
                        }
                    },
                "Chain_B":
                    {"1_GLY":
                        {"Atom_1":
                            {"atom_name":"G", 
                            "x_coordinate":3.455
                            },
                        "Atom_2":
                            {"atom_name":"G"
                            }
                        },
                    "2_SER":
                        {"Atom_1":
                            {"atom_name":"T"
                            }
                        }
                    }
                }
            }
        ]
        mycol.insert_many(mylist)