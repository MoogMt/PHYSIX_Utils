#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>
#include <typeinfo>

using namespace std;

// basic atom structure
//------------------------------
struct Atom
{
  // name
  string name;
  // position
  double x, y, z;
  // index
  int index;
};
//------------------------------


// Bond between two molecules
struct Bond
{
  int index_i, index_j;   // Index of atoms i and j
  int time;               // Nb of steps the bond lives
};

// molecules
//--------------------
struct Molecule
{
  string name;        // Name of molecule
  vector<Atom> atoms; // Atoms
};
//-------------------



// Cell
//---------------------------
struct Cell
{
  // lengths parameters
  double a, b, c;
  // angular parameters
  double alpha, beta, gamma;
};
//----------------------------


struct Contact_Matrix
{
  vector<string> types;  // Types of the atoms
  vector<double> matrix; // Contact matrix
};

// prinatoms
// Prints all informations on a given atoms in the console
//-------------------------------------------------------------------
void printAtoms( vector<Atom> atoms )
{
  for( int i=0; i < atoms.size(); i++ )
    {
      cout << "ATOM " << i+1 << " | " << atoms[i].name << " | "  << atoms[i].x << " " << " " << atoms[i].y << " " << atoms[i].z << endl;
    }
  return ;
}
//--------------------------------------------------------------------

// Min
// Returns the minimum of value of vector containing doubles
//-----------------------------------------------------------
double min(vector<double> vector)
{
  double min=vector[0];
  for( int i=0 ; i < vector.size() ; i++ )
    {
      if( min > vector[i] )
	{
	  min=vector[i];
	}
    }
  return min;
}
//-------------------------------------------

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
vector<Atom> pbc(Atom atom, Cell box)
{
  Atom atom_image;
  vector<Atom> pbc;
  
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

// Distance between two atoms
double distanceAtoms(Atom i, Atom j)
{
  double x2 = i.x - j.x;
  double y2 = i.y - j.y;
  double z2 = i.z - j.z;

  double dist=sqrt( pow(x2,2.0) + pow(y2,2.0) + pow(z2,2.0) );

  return dist;
}

// Distance
// returns the distance between two atoms in a given cell 
//------------------------------------------------------
double distanceAtoms(vector<Atom> atoms, int i, int j, Cell box)
{
  Atom atom = wrapPBC(atoms[i],box);
  vector<Atom> pbc_images = pbc(atoms[j],box);
  vector<double> distances;
  for (int l=0; l < pbc_images.size(); l++ )
    {
      distances.push_back(distanceAtoms(atom,pbc_images[l]));
    }
  return min(distances);
}
//------------------------------------------------------

// Average
// Returns the average distance between two atoms
//----------------------------------------------------------
double average(vector<int> data)
{
  double average=0;
  for( int i=0; i<data.size(); i++)   
    {
      average=average+(double)(data[i]);
    }
  return average=average/(double)(data.size()); // Returning the average
}
double average(vector<double> data)
{
  double average=0;                   // Average
  for( int i=0; i<data.size(); i++)   
    {
      average=average+data[i];
    }
  return average=average/data.size(); // Returning the average
}
//--------------------------------------

// Creates a contact matrix from a list of atoms
Contact_Matrix makeContactMatrix(vector<Atom> atom_list, Cell box)
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

int sumFromO(int integer)
{
  int sum=0;
  for ( int i=0 ; i <= integer ; i++ )
    {
      sum += i;
    }
  return sum;
}
int computeSep(int atom_index, int nb_atoms)
{
  return atom_index*nb_atoms-sumFromO(atom_index);
}

//
vector<double> getAtomContact(Contact_Matrix contact_matrix, int atom_index)
{
  vector<double> contact_atom;
  int nb_atoms = contact_matrix.types.size();
  int sep=computeSep(atom_index,nb_atoms);
  if ( atom_index != 0 )
    {
      for ( int i=0; i < atom_index; i++ )
	{
	  contact_atom.push_back(contact_matrix.matrix[atom_index+sumFromO(i)-1]);
	}
    }
  for ( int i=0; i < nb_atoms-atom_index ; i++ )
    {
      contact_atom.push_back(contact_matrix.matrix[i+sep-1]);
    }
  return contact_atom;
}

vector<int> getCoordinances(string type, Contact_Matrix contact_matrix, double cut_off_radius)
{
  vector<int> coord;
  for ( int i=0 ; i < contact_matrix.types.size() ; i++ )
    {
      if ( contact_matrix.types[i] == type )
	{
	  int neighbours=0;
	  vector<double> contact = getAtomContact(contact_matrix,i);
	  for ( int j = 0 ; j < contact.size() ; j++ )
	    {
	      if (  contact[j] < cut_off_radius)
		{
		  neighbours++;
		}
	    }
	  coord.push_back(neighbours);
	}
    }
  return coord;
}

vector<Atom> readstepXYZ(ifstream& file)
{
  // Stream Handling
  istream_iterator<string> read(file);
  istream_iterator<string> end;

  // Atoms
  int nb_atoms=-1;
  Atom atom;
  vector<Atom> atom_list;
  int atom_count=0;
  
  // Reading step
  while( read != end && atom_list.size() < nb_atoms )
    {
      if ( nb_atoms > 0 )
	{
	  // Reads one atom
	  atom.name = string(*read); ++read;
	  atom.x    = atof(string(*read).c_str()); ++read;
	  atom.y    = atof(string(*read).c_str()); ++read;
	  atom.z    = atof(string(*read).c_str());
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
	  nb_atoms = atoi(string(*read).c_str());
	  ++read; // "STEP"
	  ++read; // STEP NUMBER
	  ++read;
	}
    }

  // Return atom list
  return atom_list;
}

//================
// MAIN PROGRAM
//=====================================================================
int main(void)
{
  // Input
  //------------------------------------
  ifstream input("TRAJEC.xyz");

  //----------------------
  // Physical parameters
  //----------------------
  Cell box = {9.0,9.0,9.0,90,90,90}; // Definition of simulation box
  double cut_off_radius = 1.75;      // Cut-Off for molecules
  int step = 1;                      // Step counter
  int comp_step=5;                   // The number of step you wait to compute CM
  vector<Atom> atom_list;            // Atoms in cell
  //----------------------------

  // Reading XYZ file
  //---------------------------------------
  do
    {
      atom_list=readstepXYZ( input ); // Read one line
      if( step==2) //step % comp_step == 0 )
	{
	  Contact_Matrix contact_matrix = makeContactMatrix(atom_list,box);
	  vector<int> coordC  = getCoordinances("C",contact_matrix,cut_off_radius);
	  vector<int> coordO  = getCoordinances("O",contact_matrix,cut_off_radius);
	  cout << step << " "<< average(coordO) << endl;
	}
      step++;
    } while( atom_list.size() != 0 );
  //--------------------------------------

  input.close(); // Closing flux
 
  return 0;
}
