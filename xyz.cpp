#include "xyz.h"

std::vector<Atom> readstepXYZ(std::ifstream& file)
{
  // Stream Handling
  std::istream_iterator<std::string> read(file);
  std::istream_iterator<std::string> end;

  // Atoms
  int nb_atoms=-1;
  Atom atom;
  std::vector<Atom> atom_list;
  int atom_count=0;
  
  // Reading step
  while( read != end && atom_list.size() < nb_atoms )
    {
      if ( nb_atoms > 0 )
	{
	  // Reads one atom
	  atom.name = std::string(*read); ++read;
	  atom.x    = atof(std::string(*read).c_str()); ++read;
	  atom.y    = atof(std::string(*read).c_str()); ++read;
	  atom.z    = atof(std::string(*read).c_str());
	  atom_list.push_back(atom);
	  if ( atom_list.size() != nb_atoms )
	    {
	      ++read;
	    }
	}
      else if ( nb_atoms == 0 )
	{
	  return atom_list;
	}
      else
	{
	  nb_atoms = atoi(std::string(*read).c_str());
	  ++read; // "STEP"
	  ++read; // STEP NUMBER
	  ++read;
	}
    }

  // Return atom list
  return atom_list;
}
