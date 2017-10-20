#include "xyz.h"

//======
// READ
//==================================================================================
bool readStepXYZ( std::ifstream & file , std::vector<Atom> & atoms , std::vector<typeLUT> & lut_list , bool same_type , bool verbose )
// Reads a step of XYZ file and build atom list and LUT for types
{

  // Clearing previous values
  //---------------------------------------------------------------
  if ( atoms.size() != 0 ) atoms.clear();
  if ( same_type && lut_list.size() != 0 ) lut_list.clear();
  //---------------------------------------------------------------

  // Physical Parameters
  //---------------------------------------------------------------
  int nb_atoms = 0 ; // Numbers
  Atom atom ;        // Stock a temporary atom
  std::string line ; // Contains one line
  //---------------------------------------------------------------

  // Getting the number of atoms in the step
  //-----------------------------------------------------------------------------------
  if ( std::getline( file , line ) )
    {
      std::istringstream it_string(line);
      if ( ! ( it_string >> nb_atoms ) )
	{
	  if ( verbose ) std::cout << "Problem with file format at line 1." << std::endl;
	  return false;
	}
    }
  //-----------------------------------------------------------------------------------

  // Jumps next comment line
  //----------------------------------------------------------------------------------
  if( ! std::getline( file , line ) )
    {
      if ( verbose && line != "" )
	{
	  std::cout << "Problem with file format at line 2" << std::endl;
	}
      return false;
    }
  //----------------------------------------------------------------------------------

  // Reads atoms by atoms
  //---------------------------------------------------------------
  while ( atoms.size() < nb_atoms )
    {
      if ( std::getline( file , line ) )
	{
	  std::istringstream it_string(line);
	  if ( ! ( it_string >> atom.name >> atom.x >> atom.y >> atom.z ) )
	    {
	      if ( verbose ) std::cout << "Problem with file format at line " << 2 + atoms.size() << "."<< std::endl;
	      return false;
	    }
	  atom.index = atoms.size();
	  atoms.push_back( atom );
	  if ( ! same_type ) addAtom2LUT( lut_list , atom );
	}
    }
  //---------------------------------------------------------------
  
  return true; 
}
//-----------------------------------------------------------------------------------
std::vector<Atom> readstepXYZ(std::ifstream & file)
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
//==================================================================================

//=======
// WRITE
//==================================================================================
void writeXYZ( std::ofstream & file , std::vector<Atom> atom_list )
{
  file << atom_list.size() << std::endl;
  file << "STEP LOUTRE" << std::endl;
  for ( int i=0 ; i < atom_list.size() ; i++ )
    {
      file << atom_list[i].name << " " << atom_list[i].x <<  " " << atom_list[i].y <<  " " << atom_list[i].z << std::endl;
    }
}
//-----------------------------------------------------------------------------------
void writeXYZ( std::ofstream & file , std::vector<Atom> atom_list , int step )
{
  file << atom_list.size() << std::endl;
  file << "STEP " << step << std::endl;
  for ( int i=0 ; i < atom_list.size() ; i++ )
    {
      file << atom_list[i].name << " " << atom_list[i].x <<  " " << atom_list[i].y <<  " " << atom_list[i].z << std::endl;
    }
}
//==================================================================================
