#include "cell.h"


//
//----------------------------------------
double backIn(double x, double a)
{
  int sign;
  if( x > 0)
    {
      sign=-1;
    }
  else
    {
      sign=+1;
    }

  while( x > a || x < 0 )
    {
      x += sign*a;
    }
  return x;
}
//-----------------------------------------

Atom wrapPBC(Atom atom_in, Cell box)
{
  Atom atom_out;
  atom_out.x = backIn( atom_in.x , box.a );
  atom_out.y = backIn( atom_in.y , box.b );
  atom_out.z = backIn( atom_in.z , box.c );
  return atom_out;
}

// PBC
// Generates all the pbc image of an atom
//------------------------------------------
std::vector<Atom> pbc(Atom atom, Cell box)
{
  Atom atom_image;
  std::vector<Atom> pbc;
  
  // Wrapping PBC
  atom=wrapPBC(atom,box);
  
  pbc.push_back(atom); // Atom at the center

  // All the images
  for ( int i = -1 ; i <= 1 ; i++ )
    {
        for ( int j = -1 ; j <= 1 ; j++ )
	  {
	      for ( int k = -1 ; k <= 1 ; k++ )
		{
		  atom_image.x = atom.x + i*box.a;
		  atom_image.y = atom.y + j*box.b;
		  atom_image.z = atom.z + k*box.c;
		  pbc.push_back(atom_image);
		}
	  }
    }
  // returns the vector
  return pbc;
}

// Distance
// returns the distance between two atoms in a given cell 
//------------------------------------------------------
double distanceAtoms(std::vector<Atom> atoms, int i, int j, Cell box)
{
  Atom atom = wrapPBC(atoms[i],box);
  std::vector<Atom> pbc_images = pbc(atoms[j],box);
  std::vector<double> distances;
  for (int l=0; l < pbc_images.size(); l++ )
    {
      distances.push_back(distanceAtoms(atom,pbc_images[l]));
    }
  return min(distances);
}
//------------------------------------------------------

// READING FILES
// Custom Cell File - By Step
Cell readParamCellStep( std::ifstream& file )
{
  std::istream_iterator<std::string> read(file);
  std::istream_iterator<std::string> end;

  Cell cell;

  int count = 0;
  while( read != end && count < 6 )
    {
      switch(count)
	{
	case 0:
	  cell.a = it2real(read);
	  break;
	case 1:
	  cell.b = it2real(read);
	  break;
	case 2:
	  cell.c = it2real(read);
	  break;
	case 3:
	  cell.alpha = it2real(read);
	  break;
	case 4:
	  cell.beta  = it2real(read);
	  break;
	case 5:
	  cell.gamma = it2real(read);
	  break;
	}
      ++read;
      count++;
    }
  return cell;
}

Cell readParamCell( std::string file_name )
{

  std::ifstream file( file_name.c_str() );
  std::istream_iterator<std::string> read( file );
  std::istream_iterator<std::string> end;

  Cell cell;
  int count = 0;
  while( read != end && count < 6 )
    {
      switch(count)
	{
	case 0:
	  cell.a = it2real(read);
	  break;
	case 1:
	  cell.b = it2real(read);
	  break;
	case 2:
	  cell.c = it2real(read);
	  break;
	case 3:
	  cell.alpha = it2real(read);
	  break;
	case 4:
	  cell.beta  = it2real(read);
	  break;
	case 5:
	  cell.gamma = it2real(read);
	  break;
	}
      ++read;
      count++;
    }
  return cell;
}
