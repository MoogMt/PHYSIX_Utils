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

std::vector<int> getCoordinances(std::string type, Contact_Matrix contact_matrix, double cut_off_radius)
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

void writeCoordinance( std::ofstream& file_handle, Contact_Matrix contact_matrix,  std::string atom_type, double cut_off_radius, int step, bool alone)
{
  std::vector<int> coord_type  = getCoordinances(atom_type,contact_matrix,cut_off_radius);
  if ( !alone )
    {
      file_handle << step << " ";
    }
  for ( int i=0 ; i < coord_type.size() ; i++ )
    {
      file_handle << coord_type[0] << " ";
    }
  file_handle << average(coord_type) << " ";
  if ( !alone )
    {
      file_handle << std::endl;
    }
 
  return;
}
