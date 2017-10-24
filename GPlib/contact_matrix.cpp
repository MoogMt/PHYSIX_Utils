#include "contact_matrix.h"

//==================
// CONTACT MATRIX
//=======================================================================================
// Constructs the contact matrix
Contact_Matrix makeContactMatrix(std::vector<Atom> atom_list, Cell box)
{
  Contact_Matrix contact_matrix;
  for ( int i=0 ; i < atom_list.size() ; i++ )
    {
      contact_matrix.types.push_back(atom_list[i].name);
      for (int j=i+1; j < atom_list.size() ; j++ )
	{
	  contact_matrix.matrix.push_back( distanceAtoms( atom_list , i , j , box) );
	}
    }
  return contact_matrix;
}
//--------------------------------------------------------------------------------------
ContactMatrix makeContactMatrix ( std::vector<Atom> atom_list, std::vector<TypeLUT> lut_list , Cell box )
{
  std::vector<double> matrix;
  for ( int i=0 ; i < atom_list.size() ; i++ )
    {
      for( int j=0 ; j < atom_list.size() ; j++ )
	{
	  matrix.push_back( distanceAtoms( atom_list , i , j , box ) );
	}
    }
  return { lut_list , matrix };
}
//--------------------------------------------------------------------------------------
ContactMatrix makeContactMatrix ( AtomList atom_list, std::vector<TypeLUT> lut_list , Cell cell )
{
  std::vector<double> matrix;
  for ( int i=0 ; i < atom_list.names.size() ; i++ )
    {
      for ( int j=0 ; j < atom_list.names.size() ; j++ )
	{
	  matrix.push_back( distanceAtomsSq( atom_list , i , j , cell) );
	}
    }
  return { lut_list, matrix };
}
//--------------------------------------------------------------------------------------
ContactMatrix makeContactMatrix ( AtomList atom_list, Cell cell , CutOffMatrix cut_off , std::vector<TypeLUT> lut_list )
{
  std::vector<double> matrix;
  for ( int i=0 ; i < atom_list.names.size() ; i++ )
    {
      for ( int j=0 ; j < atom_list.names.size() ; j++ )
	{
	  int value = 0;
	  double cutoff = getCutOff( cut_off , i , j );
	  if ( distanceAtomsSq( atom_list , i , j , cell) < cutoff*cutoff )
	    {
	      value = 1;
	    }
	  matrix.push_back( value );
	}
    }
  return { lut_list, matrix };
}
//=======================================================================================

//==========
// DISTANCE
//=======================================================================================
double getDistance(Contact_Matrix contact_matrix, int atom_index_1, int atom_index_2 )
{
  int nb_atoms = contact_matrix.types.size();
  int mini=min(atom_index_1,atom_index_2);
  int maxi=max(atom_index_1,atom_index_2);
  return contact_matrix.matrix[computeSep(mini,nb_atoms)-mini+maxi-1];
}
//=======================================================================================

//==================
// ANGLES FUNCTIONS
//=======================================================================================
// Calculates the angles between three atoms centered on atom_center_index using contact matrix and Al-Kashi theorem
double getAngle( Contact_Matrix contact_matrix, int atom_center_index , int atom_2_index, int atom_3_index )
{
  double a = getDistance( contact_matrix, atom_center_index, atom_2_index);
  double b = getDistance( contact_matrix, atom_center_index, atom_3_index);
  double c = getDistance( contact_matrix, atom_center_index, atom_3_index);
  return acos( (a*a+b*b-c*c)/(2*a*b) );
}
//=======================================================================================

//================
// ATOMIC CONTACT
//=======================================================================================
std::vector<double> getAtomContact( Contact_Matrix contact_matrix , int atom_index )
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
//----------------------------------------------------------------------------------------
std::vector<double> getAtomContact( Contact_Matrix contact_matrix , int atom_index, std::string specie ) 
{
  std::vector<double> contact_atom;
  int nb_atoms = contact_matrix.types.size();
  for(int i=0 ; i < nb_atoms ; i++ )
    {
      if ( i != atom_index && contact_matrix.types[i] == specie )
	{
	  contact_atom.push_back( getDistance(contact_matrix, atom_index , i ) );
	}
    }
  return contact_atom ;
}
//-----------------------------------------------------------------------------------------
void writeAtomContact( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> atom_indexes )
// Writes the atomic contact of a set of atoms in file
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
//=======================================================================================

//=============
// COORDINANCE
//=======================================================================================
int getAtomNeighboursNb( Contact_Matrix contact_matrix, int atom_index, double cut_off_radius )
{
  int neighbours = 0;
  std::vector<double> contact = getAtomContact( contact_matrix, atom_index );
  for ( int i=0 ; i < contact.size() ; i++ )
    {
      if ( contact[i] < cut_off_radius )
	{
	  neighbours++;
	}
    }
  return neighbours;
}
//------------------------------------------------------------------------------------------------
int getAtomNeighboursNb( Contact_Matrix contact_matrix, int atom_index, std::string specie, double cut_off_radius )
{
  int neighbours = 0;
  std::vector<double> contact = getAtomContact( contact_matrix, atom_index );
  for ( int i=0 ; i < contact.size() ; i++ )
    {
      if ( contact[i] < cut_off_radius && contact_matrix.types[i] == specie )
	{
	  neighbours++;
	}
    }
  return neighbours;
}
//------------------------------------------------------------------------------------------------
std::vector<int> getAtomsNeighboursNb( Contact_Matrix contact_matrix , std::vector<int> atom_index_list , double cut_off_radius )
{
  std::vector<int> neighbours_nb_list;
  for ( int i=0 ; i < atom_index_list.size() ; i++ )
    {
      neighbours_nb_list.push_back( getAtomNeighboursNb( contact_matrix , atom_index_list[i] , cut_off_radius ) );
    }
  return neighbours_nb_list;
}
//------------------------------------------------------------------------------------------------
std::vector<int> getTypeNeighboursNb(Contact_Matrix contact_matrix, std::string type, double cut_off_radius ) 
{
  int nb_atoms = contact_matrix.types.size();
  std::vector<int> neighbours_nb_list;
  for ( int i=0 ; i <  nb_atoms ; i++ )
    {
      if ( contact_matrix.types[i] == type )
	{
	  neighbours_nb_list.push_back( getAtomNeighboursNb( contact_matrix , i , cut_off_radius) );
	}
    }
  return neighbours_nb_list;
}
//------------------------------------------------------------------------------------------------
double getTypeCoordinance( Contact_Matrix contact_matrix, std::string type, double cut_off_radius )
{
  return average( getTypeNeighboursNb( contact_matrix , type , cut_off_radius ) );
}
//=======================================================================================

//===============================
// NEAREST NEIGHBOURS FUNCTIONS
//=======================================================================================
// -> Get
//--------------------------------------------------------------------------------------------
// Gets the n_nearest neighbours distance from atom with index atom_index using a contact matrix
double getNNearest( Contact_Matrix contact_matrix , int n_nearest, int atom_index )
{
  return sortVector( getAtomContact( contact_matrix , atom_index ) , true )[n_nearest-1];
}
//------------------------------------------------------------------------------------------------
double getNNearest( Contact_Matrix contact_matrix , int n_nearest, int atom_index, std::string specie )
{
  return sortVector( getAtomContact( contact_matrix , atom_index, specie ) , true )[n_nearest-1];
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix contact_matrix , std::vector<int> n_nearest, int atom_index)
{
  std::vector<double> n_nearest_list;
  for( int i=0; i < n_nearest.size() ; i++ )
    {
      n_nearest_list.push_back( getNNearest( contact_matrix , n_nearest[i] , atom_index ) );
    }
  return n_nearest_list;
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix contact_matrix , int nearest, std::vector<int> atom_indexes)
{
  std::vector<double> n_nearest;
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      n_nearest.push_back( getNNearest( contact_matrix , nearest, atom_indexes[i] ) );
    }
  return n_nearest;
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix contact_matrix , int nearest, std::vector<int> atom_indexes , std::string specie )
{
  std::vector<double> n_nearest;
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      n_nearest.push_back( getNNearest( contact_matrix , nearest, atom_indexes[i] , specie) );
    }
  return n_nearest;
}
//------------------------------------------------------------------------------------------------
std::vector<std::vector<double> > getNNearest( Contact_Matrix contact_matrix , std::vector<int> nearest, std::vector<int> atom_indexes , std::string specie )
// Gets the n_nearest neighbours of given specie for a specific set of atoms whose with index atom_indexes
// Very ineffective, slow and not great to use
// Need to be corrected
{
  std::vector<std::vector<double> > n_nearest_table;
  for ( int i=0 ; i < nearest.size() ; i++ )
    {
      n_nearest_table.push_back( getNNearest( contact_matrix , nearest[i] , atom_indexes , specie )  );
    }
  return n_nearest_table;
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix contact_matrix , int n_nearest, std::string atom_type )
{

  std::vector<double> n_nearest_list;
  int nb_atoms = contact_matrix.types.size();
  for ( int i=0 ; i < nb_atoms ; i++ )
    {
      if( atom_type == contact_matrix.types[i] )
	{
	  n_nearest_list.push_back( getNNearest( contact_matrix , n_nearest , i ) );
	}
    }
  return n_nearest_list;
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix contact_matrix , int n_nearest, std::vector<std::string> atom_types )
{
  std::vector<double> n_nearest_list;
  for ( int i=0 ; i < atom_types.size() ; i++ )
    {
      std::vector<double> temp = getNNearest( contact_matrix , n_nearest , atom_types[i] );
      for( int j=0 ; j < temp.size() ; j++ )
	{
	  n_nearest_list.push_back(temp[i]);
	}
    }
  return n_nearest_list;
}
//------------------------------------------------------------------------------------------------
// -> Write
//------------------------------------------------------------------------------------------------
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , int n_nearest , std::vector<int> atom_indexes , int step )
{
  file << step << " ";
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      file << getNNearest( contact_matrix , n_nearest , atom_indexes[i]) << " ";
    }
  file << std::endl;
  return;
}
//------------------------------------------------------------------------------------------------
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> n_nearest, int atom_index)
// Writes the n_nearest atoms of atoms with index atom_index in file using contact matrix
{
  for( int i=0; i < n_nearest.size() ; i++ )
    {
      file << getNNearest( contact_matrix , n_nearest[i] , atom_index ) << " " ;
    }
  return;
}
//------------------------------------------------------------------------------------------------
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> nearest, std::vector<int> atom_indexes , int step)
// Gets and writes n earest neighbours distances 
{
  file << step << " ";
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      writeNearest( file , contact_matrix , nearest, atom_indexes[i] );
    }
  file << std::endl;
  return;
}
//=======================================================================================

//===========
// CONNECTED
//==================================================
bool connected( ContactMatrix cm , int i , int j )
{
  int nb_atoms = sqrt( cm.matrix.size() );
  if ( cm.matrix[ i*nb_atoms + j ] ) return true;
  else return false;
}
//=================================================