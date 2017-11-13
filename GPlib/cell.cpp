#include "cell.h"

//===========
// WRAP PBC
//====================================================================
double backIn ( double x , double a )
// Wraps a single position dimension inside PBC  
{
  // Determining sign
  int sign;
  if ( x > 0) sign = -1;
  else        sign = +1;
  // Going back inside the box
  while ( x > a || x < 0 )
    {
      x += sign*a;
    }
  // Returning
  return x;
}

//------------------------------------------------------
Atom wrapPBC(Atom atom_in, Cell box)
// Wraps a signle atom inside cell
// In:
// Out:
{ 
  Atom atom_out;
  atom_out.x = backIn( atom_in.x , box.a );
  atom_out.y = backIn( atom_in.y , box.b );
  atom_out.z = backIn( atom_in.z , box.c );
  return atom_out;
}
//--------------------------------------------------------------
std::vector<Atom> wrapPBC( std::vector<Atom> atoms , Cell cell )
// Wraps all atoms in PBC
{
  for ( int i=0 ; i < atoms.size() ; i++ )
    {
      atoms[i] = wrapPBC( atoms[i] , cell );
    }
  return atoms;
}
//--------------------------------------------------------------
void wrapPBC( AtomList & atom_list , int i , Cell cell )
// Wraps a specific atom in PBC
{
  atom_list.x[i] = backIn( atom_list.x[i], cell.a );
  atom_list.y[i] = backIn( atom_list.y[i], cell.b );
  atom_list.z[i] = backIn( atom_list.z[i], cell.c );
  return ;
}
//--------------------------------------------------------------
void wrapPBC( AtomList & atom_list , Cell cell )
// Wraps all atoms in PBC
{
  for ( int i=0 ; i < atom_list.names.size() ; i++ )
    {
      wrapPBC( atom_list , i , cell );
    }
  return;
}
//-----------------------------------------------------------------------
std::vector<Atom> pbcImages(Atom atom, Cell box)
// Generates all the pbc image of an atom
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
//===========================================================================================

//===========
// DISTANCE
//===========================================================================================
double distAtoms1D( double x1, double x2, double a )
// returns the distance between two atoms in a given cell 
{
  double dx = x1 - x2;
  if ( dx >  a*0.5 ) dx -= a;
  if ( dx < -a*0.5 ) dx += a;
  return dx*dx;
}
//------------------------------------------------------------------------------------------------
double distanceAtoms(std::vector<Atom> atoms, int i, int j, Cell box , bool wrap , bool sqrt_test )
// returns the distance ( square of the distance ) between two atoms 
{
  Atom atom_i = atoms[i] ; Atom atom_j =  atoms[j] ;
  if ( wrap )
    {
      atom_i = wrapPBC( atom_i , box ) ; 
      atom_j = wrapPBC( atom_j , box ) ;
    }
  double dist = 0 ;
  dist += distAtoms1D( atom_i.x , atom_j.x, box.a ) ;
  dist += distAtoms1D( atom_i.y , atom_j.y, box.b ) ;
  dist += distAtoms1D( atom_i.z , atom_j.z, box.c ) ;
  if ( sqrt_test ) return sqrt( dist );
  else return dist;
}
//------------------------------------------------------------------------------------------------
double distanceAtomsSq( AtomList & atom_list ,  int i , int  j , Cell cell , bool wrap )
{
  // Wrappring in PBC
  if ( wrap )
    {
      wrapPBC( atom_list , i , cell );
      wrapPBC( atom_list , j , cell );
    }
  // Calculating distance
  double dist = distAtoms1D( atom_list.x[i] , atom_list.x[j] , cell.a );
  dist += distAtoms1D( atom_list.y[i] , atom_list.y[j] , cell.b );
  dist += distAtoms1D( atom_list.z[i] , atom_list.z[j] , cell.c );
  return dist;
}
//-----------------------------------------------------------------------------------------------
void writeAtomDistances( std::ofstream & file , std::vector<Atom> atom_list , std::vector<int> atom_index, Cell box)
// Write atomic distances for all atom whose index are stored in atom_index using file pointer file.
{
  // First loop over all atoms
  for ( int i=0 ; i < atom_index.size() ; i++ )
    {
      // Loop over all other atoms
      for ( int j=0 ; j < atom_list.size() ; j++ )
 	{
	  // We don't calculate distance of an atom with itself
	  if ( atom_index[i] != j )
	    {
	      // Writing 
	      file << distanceAtoms(atom_list,atom_index[i],j,box) << std::endl;
	    }
	}
    }  
  return;
}
//===========================================================================================

//============
// MODIFY BOX
//===========================================================================================
Cell compressBox( Cell cell , double frac_a , double frac_b , double frac_c )
{
  cell.a *= frac_a;
  cell.b *= frac_b;
  cell.b *= frac_c;
  return cell;
}
//===========================================================================================

//========
// FILES
//===========================================================================================
// READ
//------------------------------------------------------------------------------------------
bool readParamCellStep( std::ifstream& file , Cell & cell  )
// Reading cell parameters from cell file, using cell pointer.
{
  // Variable
  std::string line; // line of file
  
  // Reads line...
  if ( std::getline( file, line ) )
    {
      // Checking line
      std::istringstream it_string(line);
      // Parsing line
      if ( !( it_string >> cell.a >> cell.b >> cell.c >> cell.alpha >> cell.beta >> cell.gamma) ) return false;
      else return true;
    }
  else return false;
}
//------------------------------------------------------------------------------------------
bool readParamCell( std::string file_name , Cell & cell )
// Reading cell parameters from cell file, using string
{
  //-------------------
  // Working variables
  //-------------------------------------------------
  std::ifstream file( file_name.c_str() );
  //-------------------------------------------------  

  if ( readParamCellStep( file , cell ) )
    {
      return true;
    }
  else
    {
      std::cout << "Problem reading file " << file_name << std::endl;
      return false;
    }
}
//==============================================================================================


//===================
// Reading Pressure
//==============================================================================================
bool readPressure( std::ifstream & input , double & pressure )
{
  //----------
  // Variable
  //--------------------
  pressure=0;
  std::string line;
  //---------------------

  //---------------------
  // Reading First line
  //--------------------------------------
  if ( ! std::getline( input , line ) )
    {
      return false;
    }
  //--------------------------------------
  
  //------------------
  // Reading pressure
  //--------------------------------------
  for ( int i=1 ; i <= 3 ; i++ )
    {
      if ( std::getline( input , line ) )
	{
	  std::istringstream it_string(line);
	  int j=0;
	  double stock;
	  while ( j < i )
	    {
	      stock=0;
	      if ( !( it_string >> stock ) ) return false;
	      j++;
	    }
	  pressure += stock;
	}
      else return false;
    }
  //--------------------------------------

  // Dividing by 3
  pressure /= 3 ;
  // No end of file met
  return true;
}
//==============================================================================================
