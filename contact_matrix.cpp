#include "contact_matrix.h"

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

std::vector<double> getNNearest( Contact_Matrix contact_matrix , int n_nearest, int atom_index )
{
  return sortVector(getAtomContact(contact_matrix,atom_index))[n_nearest-1];
}

void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> nearest, int atom_index)
{
  for( int i=0; i < nearest.size() ; i++ )
    {
      file << getNNearest(file,nearest[i],atom_index) << " " ;
    }
  return;
}


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

void writeNearest( std::vector<std::ofstream> files , Contact_Matrix contact_matrix , std::vector<int> nearest, std::vector<int> atom_indexes , int step)
{
  if ( files.size() == nearest.size() )
    {
      for ( int i=0 ; i < nearest.size() ; i++ )
	{
	  writeNearest( files[i] , contact_matrix , nearest[i], atom_indexes[i] );
	}
      return;
    }
  else
    {
      std::cout << "Size of the file vector and number of nearest neighbours do not match." << std::endl;
      return; 
    }
}

