#include "bonds.h"

Bond makeBond( atom_list[i] , atom_list[j] )
{
  Bond bond = { atom_list[i].index , atom_list[j].index, 1.0 };
  return bond;
}

std::vector<Bond> makeBonds(std::vector<Atom> atom_list, double cut_off_radius ) 
{
  std::vector<Bond> bonds;
  for ( int i=0 ; i < atom_list.size() ; i++ )
    {
      for ( int j=0 ; j < atom_list.size() ; j++ )
	{
	  if( distanceAtoms(  atom_list[i] , atom_list[j] ) < cut_off_radius )
	    {
	      bonds.push_back( makeBond( atom_list[i] , atom_list[j] ) );
	    }
	}
    }
  return bonds;
}

std::vector<Bond> makeBonds(Contact_Matrix contact_matrix , double cut_off_radius )
{
  std::vector<Bond> bonds;
  return bonds;
}
