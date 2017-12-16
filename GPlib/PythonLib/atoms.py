"""
Created on Thu Dec  7 03:19:35 2017

@author: moogm
"""

#=============================================================================#
class Atom:
    """
    Classe definissant la notion de projet
    Un projet est definit par:
    - son nom
    - un nombre max d'etudiants
    """
    # Constructor
    def __init__(self, name, position, type_name, type_index ):
        """Constructeur de notre classe"""
        self.name = name;
        self.position = position;
        self.type_name = type_name;
        self.type_index = type_index;
#-----------------------------------------------------------------------------#
class AtomList:
    """
    Classe definissant la notion de projet
    Un projet est definit par:
    - son nom
    - un nombre max d'etudiants
    """
    def __init__(self, names, positions, type_names, type_indexes ):
        """Constructeur de notre classe"""
        self.names = names;
        self.positions = positions;
        self.type_names = type_names;
        self.type_indexes = type_indexes;
#=======================
