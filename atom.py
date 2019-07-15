
class Atom:

    def __init__(self, atom_name, x_coordinate, y_coordinate, z_coordinate, occupancy, temperature_factor, element_symbol):
        self.atom_name = atom_name
        self.x_coordinate = x_coordinate
        self.y_coordinate = y_coordinate
        self.z_coordinate = z_coordinate
        self.occupancy = occupancy
        self.temperature_factor = temperature_factor
        self.element_symbol = element_symbol

    def get_aminoacid(self):
        pass
                       
    def __str__(self):
        return "%s - %.2f - %.2f - %.2f - %.2f - %.2f - %s" %(self.atom_name, self.x_coordinate, 
        self.y_coordinate, self.z_coordinate, self.occupancy, self.temperature_factor,
        self.element_symbol)
"""

    def general_dictionary(self):
        protein_dictionary = {"_id":self.protein_name}
        value_dictionary = 0

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

                    protein_dictionary[chain_identifier + "_" + residue_name + "_" + atom_serial_number]= {
                        chain_identifier:{
                            residue_name:{
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
                                }, 
                            "residue_sequence_number": residue_sequence_number, 
                            "residue_name":residue_name
                        }
                    }

            return protein_dictionary

dictionary = Protein("2ki5").general_dictionary()
collection.insert_one(dictionary)
print(dictionary)
"""