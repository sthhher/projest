
class Chain:
    
    def __init__(self, chain_identifier, protein_name):
        self.protein_name = protein_name
        self.chain_identifier = chain_identifier

    def get_chains(self, protein_name):
        with open(protein_name, 'r+') as text_file:
            for line in text_file:
                if line.startswith("ATOM"):
                    self.chain_identifier = line[21]
                    self.chain_identifier = self.chain_identifier.strip()

                    return self.chain_identifier

    def __str__(self):
        return "%s" %(self.chain_identifier)