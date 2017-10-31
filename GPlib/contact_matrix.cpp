#include "contact_matrix.h"

//==================
// CONTACT MATRIX
//=======================================================================================
Contact_Matrix makeContactMatrix(std::vector<Atom> atom_list, Cell box)
// Constructs the restricted contact matrix
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
ContactMatrix makeContactMatrix ( AtomList atom_list, Cell cell , CutOffMatrix cut_off , AllTypeLUT lut_type )
// Constructs the full contact matrix
{
  // Determines the number of atoms
  int nb_atoms = atom_list.x.size();

  // Initialize contact matrix with 0s
  std::vector<double> matrix; matrix.assign(nb_atoms*nb_atoms,0.);

  // Loop over all pairs of atoms
  for ( int i=0 ; i < nb_atoms-1 ; i++ )
    {
      for ( int j=i+1 ; j < nb_atoms ; j++ )
	{
	  // Getting the cut_off for the atom_i vs atom_j interaction
	  double cutoff = getCutOff( cut_off , lut_type.type_index[i] , lut_type.type_index[j] );
	  // Comparing distance to cut_off
	  if ( distanceAtomsSq( atom_list , i , j , cell) < cutoff*cutoff )
	    {
	      matrix[i*nb_atoms+j] = 1;
	      matrix[j*nb_atoms+i] = 1; 
	    }
	}
    }

  // Sending results
  return { nb_atoms, lut_type, matrix };
}
//=======================================================================================

//===========
// CONNECTED
//==================================================
bool connected( ContactMatrix cm , int i , int j )
{
  if ( cm.matrix[ i*cm.nb_atoms + j ] ) return true;
  else return false;
}
//=================================================

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
//--------------------------------------------------------------------------------------
double getDistance( ContactMatrix cm, int atom_index1 , int atom_index2 )
{
  return cm.matrix[atom_index1*cm.nb_atoms+atom_index2];
}
//=======================================================================================

//==================
// ANGLES FUNCTIONS
//=======================================================================================
// Calculates the angles between three atoms centered on atom_center_index using contact matrix and Al-Kashi theorem
double getAngle( Contact_Matrix contact_matrix, int atom_center_index , int atom_2_index, int atom_3_index )
// Get Angle
{
  double a = getDistance( contact_matrix, atom_center_index, atom_2_index);
  double b = getDistance( contact_matrix, atom_center_index, atom_3_index);
  double c = getDistance( contact_matrix, atom_center_index, atom_3_index);
  return acos( (a*a+b*b-c*c)/(2*a*b) );
}
//--------------------------------------------------------------------------------------
double getAngle( ContactMatrix cm , int atom_A , int atom_B , int atom_C )
// Get Angle(A)=(BAC)
{
  double a = getDistance( cm , atom_A , atom_C );
  double b = getDistance( cm , atom_A , atom_C );
  double c = getDistance( cm , atom_B , atom_C );
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
std::vector<double> getAtomContact( ContactMatrix cm , int atom_index )
{
  std::vector<double> contact;
  int offset = cm.nb_atoms*atom_index;
  for ( int i=0 ; i < cm.nb_atoms ; i++ )
    {
      contact.push_back(cm.matrix[offset+i]);
    }
  return contact;
}
//-----------------------------------------------------------------------------------------
std::vector<double> getAtomContact( ContactMatrix cm , int atom_index , std::string specie )
{
  std::vector<double> contact;
  int offset = cm.nb_atoms*atom_index;
  for ( int i=0 ; i < cm.nb_atoms ; i++ )
    {
      if ( cm.lut_list.type_name[atom_index] == specie )
	{
	  contact.push_back(cm.matrix[offset+i]);
	}
    }
  return contact;
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
// RESTRICTED CONTACT MATRIX
//------------------------------------------------------------------------------------------------
double getNNearest( Contact_Matrix cm , int n_nearest, int atom_index )
// Gets the n_nearest neighbours distance from atom with index atom_index using a contact matrix
{
  return sortVector( getAtomContact( cm , atom_index ) , true )[n_nearest-1];
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix cm , std::vector<int> n_nearest, int atom_index )
{
  std::vector<double> n_nearest_list;
  for( int i=0; i < n_nearest.size() ; i++ )
    {
      n_nearest_list.push_back( getNNearest( cm , n_nearest[i] , atom_index ) );
    }
  return n_nearest_list;
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix cm , int nearest, std::vector<int> atom_indexes )
{
  std::vector<double> n_nearest;
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      n_nearest.push_back( getNNearest( cm , nearest, atom_indexes[i] ) );
    }
  return n_nearest;
}
//------------------------------------------------------------------------------------------------
// By Specie
//------------------------------------------------------------------------------------------------
double getNNearest( Contact_Matrix cm , int n_nearest, int atom_index, std::string specie )
{
  return sortVector( getAtomContact( cm , atom_index, specie ) , true )[n_nearest-1];
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix cm , int n_nearest , std::string specie )
{
  return sortVector( getAtomContact( cm , n_nearest , specie ) , true );
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix cm , int n_nearest, std::vector<std::string> atom_types )
{
  std::vector<double> n_nearest_list;
  for ( int i=0 ; i < atom_types.size() ; i++ )
    {
      std::vector<double> temp = getNNearest( cm , n_nearest , atom_types[i] );
      for( int j=0 ; j < temp.size() ; j++ )
	{
	  n_nearest_list.push_back(temp[i]);
	}
    }
  return n_nearest_list;
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix cm , int nearest, std::vector<int> atom_indexes , std::string specie )
{
  std::vector<double> n_nearest;
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      n_nearest.push_back( getNNearest( cm , nearest, atom_indexes[i] , specie) );
    }
  return n_nearest;
}
//------------------------------------------------------------------------------------------------
// FULL CONTACT MATRIX
//------------------------------------------------------------------------------------------------
double getNNearest( ContactMatrix cm , int n_nearest , int atom_index )
{
  return sortVector( getAtomContact( cm , atom_index ) , true )[n_nearest-1];
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( ContactMatrix cm , std::vector<int> n_nearest, int atom_index)
{
  std::vector<double> n_nearest_list;
  for( int i=0; i < n_nearest.size() ; i++ )
    {
      n_nearest_list.push_back( getNNearest( cm, n_nearest[i] , atom_index ) );
    }
  return n_nearest_list;
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( ContactMatrix cm , int nearest, std::vector<int> atom_indexes)
{
  std::vector<double> n_nearest;
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      n_nearest.push_back( getNNearest( cm , nearest , atom_indexes[i] ) );
    }
  return n_nearest;
}
//------------------------------------------------------------------------------------------------
// By Specie
//------------------------------------------------------------------------------------------------
double getNNearest( ContactMatrix cm , int atom_index , int n_nearest , std::string specie )
{
  return sortVector( getAtomContact( cm , atom_index, specie ) , true )[n_nearest-1];
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearestVector( ContactMatrix cm , int atom_index , int n_nearest , std::string specie )
{
  return sortVector( getAtomContact( cm , atom_index, specie ) , true );
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( ContactMatrix cm , int nearest, std::vector<int> atom_indexes , std::string specie )
{
  std::vector<double> n_nearest;
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      n_nearest.push_back( getNNearest( cm , nearest, atom_indexes[i] , specie) );
    }
  return n_nearest;
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( ContactMatrix cm , int atom_index, int n_nearest, std::vector<std::string> atom_types )
{
  std::vector<double> n_nearest_list;
  for ( int i=0 ; i < atom_types.size() ; i++ )
    {
      std::vector<double> temp = getNNearestVector( cm , atom_index , n_nearest , atom_types[i] );
      for( int j=0 ; j < temp.size() ; j++ )
	{
	  n_nearest_list.push_back(temp[i]);
	}
    }
  return n_nearest_list;
}
//===================================================================================================

//======
// IO
//==================================================================================================
// PRINT
//------------------------------------------------------------------------------------------------
void printContactMatrix( ContactMatrix cm )
{
  for (int i=0; i < cm.nb_atoms*cm.nb_atoms; i++)
    {
      std::cout << cm.matrix[i] << " ";
      if ( (i+1) % cm.nb_atoms == 0 )
	{
	  std::cout << std::endl;
	}
    }
}
//------------------------------------------------------------------------------------------------
// WRITE
//------------------------------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------------------
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
//==================================================================================================
