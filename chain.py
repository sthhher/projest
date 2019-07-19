
class Chain:
    
    def __init__(self, chain_identifier, protein_name):
        self.protein_name = protein_name
        self.chain_identifier = chain_identifier
        
    def __str__(self):
        return "%s" %(self.chain_identifier)