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
  double x;
  double y;
  double z;
  // index
  int index;
};
//------------------------------

// molecules
//--------------------
struct molecule
{
  // atom list
  vector<Atom> list;
};
//-------------------

// Cell
//---------------------------
struct cell
{
  // lengths parameters
  double a, b, c;
  // angular parameters
  double alpha, beta, gamma;
};
//----------------------------


struct Contact_Matrix
{
  vector<string> types;  // types of the atoms
  vector<double> matrix; // contact matrix
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

Atom wrapPBC(Atom atom_in, cell box)
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
vector<Atom> pbc(Atom atom, cell box)
{
  Atom atom_image;
  vector<Atom> pbc;
  
  // Wrapping PBC
  atom=wrapPBC(atom,box);
  
  pbc.push_back(atom); // Atom at the center

  // All the images
  for ( int i = -1 ; i <= 1 ; i=i+2 )
    {
        for ( int j = -1 ; j <= 1 ; j=j+2 )
	  {
	      for ( int k = -1 ; k <= 1 ; k=k+2 )
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
  double x = i.x - j.x;
  double y = i.y - j.y;
  double z = i.z - j.z;

  double dist=sqrt( pow(x,2.0) + pow(y,2.0) + pow(z,2.0) );
  if ( dist > 9 )
    {
      cout << "AHAH" << endl;
      cout << "x = " << x << endl;
      cout << "x = "<< i.x << " x'=" << j.x << endl;
      cout << "y = " << y << endl;
      cout << "y = "<< i.y << " y'=" << j.y << endl;
      cout << "z = " << z << endl;
      cout << "z = "<< i.z << " z'=" << j.z << endl;
    }
  return dist;
}

// Distance
// returns the distance between two atoms in a given cell 
//------------------------------------------------------
double distanceAtoms(vector<Atom> atoms, int i, int j, cell box)
{
  vector<Atom> pbc_images = pbc(atoms[j],box);
  vector<double> distances;
  for (int i=0; i < pbc_images.size(); i++ )
    {
      distances.push_back(distanceAtoms(wrapPBC(atoms[i],box),pbc_images[i]));
    }
  return min(distances);
  }
//------------------------------------------------------

// Average
// Returns the average distance between two atoms
//----------------------------------------------------------
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
Contact_Matrix makeContactMatrix(vector<Atom> atom_list, cell box)
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

//
int computeSep(int atom_index,int nb_atoms)
{
  int sep=0;
  for (int i=1 ; i < atom_index+1; i++ )
    {
      sep += nb_atoms -i;
    }
  return sep;
}

//
vector<double> getAtomContact(Contact_Matrix contact_matrix, int atom_index)
{
  vector<double> contact_atom;
  int nb_atoms = contact_matrix.types.size();
  int sep = computeSep(atom_index,nb_atoms);
  for ( int i=1; i < sep; i=i+nb_atoms-2 )
    {
      contact_atom.push_back(contact_matrix.matrix[i]);
    }
  for ( int i=sep; i < sep + nb_atoms - atom_index ; i++ )
    {
      contact_atom.push_back(contact_matrix.matrix[i]);
    }
  return contact_atom;
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
  cell box = {9.0,9.0,9.0,90,90,90};
  double cut_off_radius = 1.75;
  int n_atoms =  96;
  int step=1;
  vector<Atom> atom_list;
  //----------------------------

  // Reading XYZ file
  //---------------------------------------
  do
    {
      atom_list=readstepXYZ( input );
      if( step % 5 == 0 )
	{
	  Contact_Matrix contact_matrix=makeContactMatrix(atom_list,box);
	  vector<double> contact=getAtomContact(contact_matrix,1);
	  for ( int i=contact.size()-64 ; i < contact.size() ; i++ )
	    {
	      if ( contact[i] > 5 )
		{
		  //		  cout << step << " " << contact[i] << " " << i << endl;
		}
	    }
	  cout << endl;
	}
      step++;
    } while( atom_list.size() != 0 );
  //--------------------------------------

  // Closing flux
  input.close();
 
  return 0;
}
