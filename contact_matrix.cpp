#include "contact_matrix.h"

//----------------
// CONTACT MATRIX
//---------------------------------------------------------------------------------
// Constructs the contact matrix
Contact_Matrix makeContactMatrix(std::vector<Atom> atom_list, Cell box)
{
  Contact_Matrix contact_matrix;
  for ( int i=0 ; i < atom_list.size() ; i++ )
    {
      contact_matrix.types.push_back(atom_list[i].name);
      for (int j=i+1; j < atom_list.size() ; j++ )
	{
	  contact_matrix.matrix.push_back(distanceAtoms(atom_list,i,j,box));
	}
    }
  return contact_matrix;
}
//---------------------------------------------------------------------------------

//-----------
// DISTANCES
//----------------------------
double getDistance(Contact_Matrix contact_matrix, int atom_index_1, int atom_index_2 )
{
  int nb_atoms = contact_matrix.types.size();
  int mini=min(atom_index_1,atom_index_2);
  int maxi=max(atom_index_1,atom_index_2);
  return contact_matrix.matrix[computeSep(mini,nb_atoms)+maxi];
}
//----------------------------

//-----------------
// ANGLES FUNCTIONS
//-------------------------------------------------------------------------------------------------------------------
// Calculates the angles between three atoms centered on atom_center_index using contact matrix and Al-Kashi theorem
double getAngle( Contact_Matrix contact_matrix, int atom_center_index , int atom_2_index, int atom_3_index )
{
  double a = getDistance( contact_matrix, atom_center_index, atom_2_index);
  double b = getDistance( contact_matrix, atom_center_index, atom_3_index);
  double c = getDistance( contact_matrix, atom_center_index, atom_3_index);
  return acos( (a*a+b*b-c*c)/(2*a*b) );
}
//-------------------------------------------------------------------------------------------------------------------

//-----------------
// ATOMIC CONTACT
//--------------------------------------------------------------------------------------------------
std::vector<double> getAtomContact(Contact_Matrix contact_matrix, int atom_index)
{
  std::vector<double> contact_atom;
  int nb_atoms = contact_matrix.types.size();
  int sep=computeSep(atom_index,nb_atoms);
  for ( int i=0 ; i < atom_index ; i++ )
    {
      contact_atom.push_back(contact_matrix.matrix[atom_index+sumBtw(nb_atoms-1,nb_atoms-i-1)-1]);
    }
  for ( int i=0; i < nb_atoms-atom_index-1 ; i++ )
    {
      contact_atom.push_back(contact_matrix.matrix[sep+i]);
    }
  return contact_atom;
}
// Writes the atomic contact of a set of atoms in file
void writeAtomContact( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> atom_indexes )
{

  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      std::vector<double> contact = getAtomContact(contact_matrix,atom_indexes[i]);
      
      for ( int j=0 ; j < contact.size() ; j++ )
	{
	  file << contact[j] << " " ;
	}
    }
  file << std::endl;
  return;
}
//---------------------------------------------------------------------------------------------------

//---------------
// COORDINANCE
//-----------------------------------------------------------------------------------------------------------
std::vector<int> getTypeCoordinances(std::string type, Contact_Matrix contact_matrix, double cut_off_radius)
{
  std::vector<int> coord;
  for ( int i=0 ; i < contact_matrix.types.size() ; i++ )
    {
      if ( contact_matrix.types[i] == type )
	{
	  int neighbours=0;
	  std::vector<double> contact = getAtomContact( contact_matrix, i );
	  for ( int j = 0 ; j < contact.size() ; j++ )
	    {
	      if ( contact[j] < cut_off_radius )
		{
		  neighbours++;
		}
	    }
	  coord.push_back(neighbours);
	}
    }

  return coord;
}
//----------------------------------------------------------------------------------------------------------------

//--------------------------------
// NEAREST NEIGHBOURS FUNCTIONS
//-------------------------------------------------------------------------------------------------
// Gets the n_nearest neighbours distance from atom with index atom_index using a contact matrix
double getNNearest( Contact_Matrix contact_matrix , int n_nearest, int atom_index )
{
  return sortVector(getAtomContact(contact_matrix,atom_index),true)[n_nearest-1];
}
// Writes the n_nearest atoms of atoms with index atom_index in file using contact matrix
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> n_nearest, int atom_index)
{
  for( int i=0; i < n_nearest.size() ; i++ )
    {
      file << getNNearest( contact_matrix , n_nearest[i] , atom_index ) << " " ;
    }
  return;
}
// Gets and writes n earest neighbours distances 
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> nearest, std::vector<int> atom_indexes , int step)
{
  file << step << " ";
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      writeNearest( file , contact_matrix , nearest, atom_indexes[i] );
    }
  file << std::endl;
  return;
}
//----------------------------------------------------------------------------------------

//----------------------
// OTHER (TO BE MOVED)
//-------------------------------------------------------------------------------------------------------------------
void writeAtomDistances( std::ofstream & file , std::vector<Atom> atom_list , std::vector<int> atom_index, Cell box)
{
  for ( int i=0 ; i < atom_index.size() ; i++ )
    {
      for ( int j=0 ; j < atom_list.size() ; j++ )
 	{
	  if ( atom_index[i] != j )
	    {
	      file << distanceAtoms(atom_list,atom_index[i],j,box) << std::endl;
	    }
	}
    }
  file << std::endl;
  return;
}
//-------------------------------------------------------------------------------------------------------------------
