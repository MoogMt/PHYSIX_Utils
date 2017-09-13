#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

using namespace std;


// basic atom structure
//------------------------------
struct atom
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
  vector<atom> list;
};
//-------------------

// Cell
//---------------------------
struct cell
{
  double a, b, c;
  double alpha, beta, gamma;
};
f//----------------------------

// prinatoms
// Prints all informations on a given atoms in the console
//-------------------------------------------------------------------
void printatoms( vector<atom> atoms )
{
  cout << "--------------------------------------------" << endl;
  for( int i=0; i < atoms.size(); i++ )
    {
      cout << "ATOM " << i+1 << " | " << atoms[i].name << " | "  << atoms[i].x << " " << " " << atoms[i].y << " " << atoms[i].z << endl;
    }
  return ;
  cout << "-------------------------------------------" << endl;
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

// PBC
// returns the shortest squared distance in 1D taking acount of pbc a
//------------------------------------------
double pbc(double x1, double x2, double a)
{
  vector<double> x;
  x.push_back( pow( x1-x2   ,2) );
  x.push_back( pow( x1+a-x2 ,2) );
  x.push_back( pow( x1-a-x2 ,2) );
  return min(x);
}

// Distance
// returns the distance between two atoms in a given cell 
//------------------------------------------------------
double distance_atoms(vector<atom> atoms, int i, int j, double a, double b, double c)
{
  double x=pbc(atoms[i].x,atoms[j].x,a);
  double y=pbc(atoms[i].y,atoms[j].y,b);
  double z=pbc(atoms[i].z,atoms[j].z,c);
  return sqrt(x+y+z);
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

// Make molecules
// returns a molecule group using a list of atom and a given cut-off radius
//-------------------------------------------------------
vector<molecule> make_molecules( vector<atom> atom_list, double cut_off_radius)
{
  vector<molecule> molecule_group; // Storing the molecules
  molecule molecule_temp;          // Temp storage for molecules
  // Loop over atoms
  for ( int i=0 ; i < atom_list.size() ; i++ )
    {
      if ( atom_list[i].name.compare("C") == 0 )
	{
	  molecule_temp.list.push_back(atom_list[i]); // Adding the C to the molecule
	  for ( int j=0 ; j < atom_list.size() ; j++)
	    {
	      if ( (atom_list[j].name).compare("O") == 0 && distance_atoms(atom_list,i,j,10.112,10.112,10.112) < cut_off_radius && i!=j )
		{
		  molecule_buff.list.push_back(atom_list[j]); // Adding the O atoms
		}
	    }
	  molecules_group.push_back(molecule_buff);
	  molecule_temp.list.clear(); // Clearing the buffer
	}
    }
  return molecule_group; // Return the molecules group
}
//-----------------------------------------------------


// GetMoleculesWithN
// Extract the molecules that have only N atoms from a molecule group
//--------------------------------------------------------
vector<molecule> getMoleculeWithN(vector<molecule> molecule_group, int n)
{
  vector<molecule> molecule_n_atoms;            // Group of molecules with only n atoms
  for( int i=0; i<molecule_group.size(); i++)   // For each molecules in the group
    {
      if( (molecule_group[i].list).size() == n) // If the molecule has n atoms...
	{
	  molecule_n_atoms.push_back(mols[i]); // ... add molecule to the group
	}
    }
  return mols2; // ... and finally return the group
}
//--------------------------------------------------------

// distance_avg_mol
// Return the average distance for a given molecule
//--------------------------------------------------------------------------
double distance_avg_mol(molecule mol, double a, double b, double c)
{
  double distance=0;
  for ( int i=1; i < (mol.list).size() ; i ++ )
    {
      distance =+ distance2(mol.list,0,i,a,b,c);
    }
  return distance/(double)(mol.list.size());
}
//--------------------------------------------------------------------------

// 
// Returns the average distances for a group of molecules
//-----------------------------------------------------------------------------
double distance_avg_mols(vector<molecule> molecule_group,double a, double b, double c)
{
  double distance=0;
  int count = 0;
  for( int i=0; i < mols.size() ; i++)
    {
      for ( int j=1; j < mols[i].list.size(); j++ )
	{
	  distance =+ distance2(mols[i].list,0,j,a,b,c);
	}
    }
  return distance/(double)count;
}
//----------------------------------------------------------------------------


// Gives the average distance for a molecule of with a given number of atoms
//---------------------------------------------------------------------------------------
double distance_avg_molsn(vector<molecule> molecules_group, int n, double a, double b, double c)
{
  vector<molecules> molecules_n_atoms = getMoleculeWithN(vector<molecule> molecule_group, int n);
  if( molecules_n_atoms.size() != 0 )
    {
      return distance_avg_molecules(molecules_n_atoms,a,b,c)
    }
  else
    {
      return 0;
    }
}
//---------------------------------------------------------------------------------------

//================
// MAIN PROGRAM
//=====================================================================
int main(void)
{
  // Atoms and molecules
  //-----------------------
  atom atom; // Temp variable that creates atoms
  vector<atom> atom_list; // Temp Variable that contains the atoms of the system
  vector<molecule> mols;  // Contains all the molecules of the system
  vector<molecule> mols_mem; // Contains all the molecules of the system at step 1
  vector<molecule> mols2;    // Contains all molecules with 2 atoms
  vector<molecule> mols3;    //  ""      ""     ""       "" 3 "" 
  vector<molecule> mols4;    //  ""      ""     ""       "" 4 "" 
  //-----------------------
  
  // Input
  //------------------------------------
  ifstream fichier("TRAJEC.xyz");
  // Input Flux Iterator
  istream_iterator<string> it(fichier);
  istream_iterator<string> end;
  //------------------------------------

  // Output
  //-------------
  ofstream length;
  ofstream COn;
  ofstream diff;
  ofstream diffper;
  // Opening files
  //--------------------------
  length.open ("length.dat");
  COn.open ("COn.dat");
  diff.open("diff.dat");
  diffper.open("diff_per.dat");
  //--------------------------

  //---------------

  // CELL length
  double a=9.0, b=9.0, c=9.0;
  // Stuff
  double cut_off_radius = 1.75;
  
   
  // Count
  //---------------------------------------------
  int n=0;           // Number of elements
  int step=0;        // Step number
  int poscount=0;    // Position count for x,y,z
  int atom_count=1;  // number of atoms
  //----------------------------------------------

  // Data
  //vector<double> angles;
  double diff_dist;
  //
   
  // Reading the file 
  while(it != end )   
    {
      //------------------------
      // Item number iteration
      //------------------------
      ++n;
      //------------------------
      if ( n > 3 )
	{
	  if( string(*it).compare("C") == 0 ||  string(*it).compare("O") == 0 )
	    {
	      atom.name=string(*it);
	    }
	  else
	    {
	      poscount++;
	      switch(poscount)
		{
		case 1:
		  atom.index=atom_count;
		  atom.x = atof(string(*it).c_str());
		  break;
		case 2:
		  atom.y = atof(string(*it).c_str());
		  break;
		case 3:
		  atom.z = atof(string(*it).c_str());
		  atom_count++;
		  atom.push_back(atom_buff);
		  poscount=0;
		  break;
		default:
		  cout << "Something went terribly wrong" << endl;
		  return 0;
		  break;	  
		}
	    }
	}
      //-----------------
      // Compiling data
      //--------------------------------------------------------
      int count_bond=0;
      int count_diff=0;
      if ( n == 32*3*4+3 )
	{
	  mols=make_molecule(atom_list); // Making molecules from atom list
	  if( step == 1 ) // Bookeeping original molecules for diffusion 
	    {
	      mols_mem=mols;
	    }
	  diff_dist  = 0; // Reinitialize diffusion distance
	  count_bond = 0; // Reinitialize the number of bonds in the system
	  for( int i=0; i < mols_mem.size() ; i++)
	    {
	      for( int j=1; j < mols_mem[i].list.size(); j++)
		{
		  double d=distance_atoms(atom_list,mols_mem[i].list[0].index-1,mols_mem[i].list[j].index-1,a,b,c);
		  if ( d > 1.75)
		    {
		      count_diff++;
		    }
		  diff_dist =+ d;
		  count_bond++;
		}
	    }
	  // 
	  diff_dist=diff_dist/((double)count_bond);
	  diffper << step << " " << count_diff << endl;
	  diff << step << " " << diff_dist << endl;
	  // Get different types of molecules
	  mols2=getCOn(mols,3);
	  mols3=getCOn(mols,4);
	  mols4=getCOn(mols,5);
	  // Output coordinance values
	  COn << step << " " << mols.size() << " " << mols2.size() << " " << mols3.size() << " " << mols4.size() << endl;
	  // Distances count for different molecules
	  length << step << " " <<  distance_avg_mols(mols2,2) << " " << distance_avg_mols(mols3,3) << " " <<  distance_avg_mols(mols) << " " << distance_avg_mols(mols4,4)  << endl;
	  // Iteration stuff
	  atom_count=1; // Reinitialize count for atoms
	  atom_list.clear(); // Clear the list of atoms, for next step
	  step++;  // Incrementing step count
	  n=0;     // Reinitialize the element number for reading
	}
      ++it; // Reads next element
    }
  // Closing fluxes
  length.close();
  COn.close();
  diff.close();
  diffper.close();
  // END 
  return 0;
}
