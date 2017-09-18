#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

# define BUFF_ELEMENTS 3
# define NB_ELEMENTS_LINE 4

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
  // atom list
  vector<atom> list;
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
double distance_atoms(vector<atom> atoms, int i, int j, cell box)
{
  double x=pbc(atoms[i].x,atoms[j].x,box.a);
  double y=pbc(atoms[i].y,atoms[j].y,box.b);
  double z=pbc(atoms[i].z,atoms[j].z,box.c);
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

//
Contact_Matrix make_contact_matrix(vector<atom> atom_list, cell box)
{
  Contact_Matrix contact_matrix;
  for ( int i=0 ; i < atom_list.size() ; i++ )
    {
      contact_matrix.types.push_back(atom_list[i].name);
      for (int j=i+1; i < atom_list.size() ; j++ )
	{
	  contact_matrix.matrix.push_back(distance_atoms(atom_list,i,j,box));
	}
    }
  return contact_matrix;
}

// Make molecules
// returns a molecule group using a list of atom and a given cut-off radius
//-------------------------------------------------------
vector<molecule> make_molecules( vector<atom> atom_list, double cut_off_radius, cell box)
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
	      if ( (atom_list[j].name).compare("O") == 0 && distance_atoms(atom_list,i,j,box) < cut_off_radius && i!=j )
		{
		  molecule_temp.list.push_back(atom_list[j]); // Adding the O atoms
		}
	    }
	  molecule_group.push_back(molecule_temp);
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
  vector<molecule> molecules_n_atoms;            // Group of molecules with only n atoms
  for( int i=0; i<molecule_group.size(); i++)   // For each molecules in the group
    {
      if( (molecule_group[i].list).size() == n) // If the molecule has n atoms...
	{
	  molecules_n_atoms.push_back(molecule_group[i]); // ... add molecule to the group
	}
    }
  return molecules_n_atoms; // ... and finally return the group
}
//--------------------------------------------------------

// distance_avg_mol
// Return the average distance for a given molecule
//--------------------------------------------------------------------------
double distance_avg_mol(molecule mol, cell box)
{
  double distance=0;
  for ( int i=1; i < (mol.list).size() ; i ++ )
    {
      distance += distance_atoms(mol.list,0,i,box);
    }
  return distance/(double)(mol.list.size());
}
//--------------------------------------------------------------------------

// 
// Returns the average distances for a group of molecules
//-----------------------------------------------------------------------------
double distance_avg_mols(vector<molecule> molecule_group,cell box)
{
  double distance=0;
  int count = 0;
  for( int i=0; i < molecule_group.size() ; i++)
    {
      for ( int j=1; j < molecule_group[i].list.size(); j++ )
	{
	  distance += distance_atoms(molecule_group[i].list,0,j,box);
	  count++;
	}
    }
  return distance/(double)count;
}
//----------------------------------------------------------------------------


// Gives the average distance for a molecule of with a given number of atoms
//---------------------------------------------------------------------------------------
double distance_avg_molsn(vector<molecule> molecules_group, int n, cell box)
{
  vector<molecule> molecules_n_atoms = getMoleculeWithN(molecules_group,n);
  if( molecules_n_atoms.size() != 0 )
    {
      return distance_avg_mols(molecules_n_atoms,box);
    }
  else
    {
      return 0;
    }
}
//---------------------------------------------------------------------------------------

// Returns the nearest neighbour distance
//-------------------------------------------------
double getNN_distance(vector<atom> atom_list, int index, cell box)
{
  int start =0;
  if( index == 0 )
    {
      start=1;
    }
  double distance = distance_atoms(atom_list,index,start,box);
  for( int i=start+1 ; i < atom_list.size() ; i++ )
    {
      if( i != index )
	{
	  if( distance > distance_atoms(atom_list,index,i,box) )
	    {
	      distance = distance_atoms(atom_list,index,i,box);
	    }
	}
    }
  return distance;
}
//--------------------------------------------------

double getNN2_distance(vector<atom> atom_list, int index, cell box)
{
  double distance_NN=getNN_distance(atom_list,index,box);
  int start =0;
  if( index == 0 )
    {
      start=1;
    }
  double distance = distance_atoms(atom_list,index,start,box);
  for( int i=start+1 ; i < atom_list.size() ; i++ )
    {
      if( i != index )
	{
	  if( distance > distance_atoms(atom_list,index,i,box) && distance_atoms(atom_list,index,i,box) > distance_NN )
	    {
	      distance = distance_atoms(atom_list,index,i,box);
	    }
	}
    }
  return distance;
}

//================
// MAIN PROGRAM
//=====================================================================
int main(void)
{
  // Atoms and molecules
  //-----------------------
  atom atom_temp ; // Temp variable that creates atoms
  vector<atom> atom_list ; // Temp Variable that contains the atoms of the system
  vector<molecule> molecules ;  // Contains all the molecules of the system
  vector<molecule> molecules_legacy ; // Contains all the molecules of the system at step 1
  vector<molecule> molecules1 ;    // Contains all molecules with 2 atoms
  vector<molecule> molecules2 ;    // Contains all molecules with 2 atoms
  vector<molecule> molecules3 ;    //  ""      ""     ""       "" 3 "" 
  vector<molecule> molecules4 ;    //  ""      ""     ""       "" 4 ""
  vector<molecule> molecules5 ;    // Contains all molecules with 2 atoms
  //-----------------------
  
  // Input
  //------------------------------------
  ifstream infile("TRAJEC.xyz");
  // Input Flux Iterator
  istream_iterator<string> it(infile);
  istream_iterator<string> end;
  //------------------------------------

  // Output
  //----------------------------------------------------
  ofstream avg_dist;   avg_dist.open ("length.dat");
  ofstream coordC;     coordC.open ("coordC.dat");
  ofstream coordO;     coordO.open ("coordO.dat");
  ofstream diff;       diff.open("diff.dat");
  ofstream firstNN_C;  firstNN_C.open("firstNN_C.dat");
  ofstream firstNN_O;  firstNN_O.open("firstNN_O.dat");
  ofstream secondNN_C;  secondNN_C.open("secondNN_C.dat");
  ofstream secondNN_O;  secondNN_O.open("secondNN_O.dat");
  //----------------------------------------------------

  //----------------------
  // Physical parameters
  //----------------------
  // cell parameters
  cell box = {9.0,9.0,9.0,90,90,90};
  // Stuff
  double cut_off_radius = 1.75;
  int n_atoms =  96;
  //----------------------------
  
  // Count
  //---------------------------------------------
  int n=0;           // Number of elements
  int step=0;        // Step number
  int poscount=0;    // Position count for x,y,z
  int atom_count=0;  // number of atoms
  int count_bond=0;  // number of bonds
  int count_diff=0;  // number of diffused O
  //----------------------------------------------

  // Data
  //vector<double> angles;
  double diff_dist;
   
  // Reading the file 
  while(it != end )   
    {
      // Item number iteration
      ++n;
      if ( n > BUFF_ELEMENTS )
	{
	  if( string(*it).compare("C") == 0 ||  string(*it).compare("O") == 0 )
	    {
	      atom_temp.name=string(*it);
	    }
	  else
	    {
	      poscount++;
	      switch(poscount)
		{
		case 1:
		  atom_temp.index=atom_count; // Atom number 
		  atom_temp.x = atof(string(*it).c_str()); // atom x position
		  break;
		case 2:
		  atom_temp.y = atof(string(*it).c_str()); // atom y position
		  break;
		case 3:
		  atom_temp.z = atof(string(*it).c_str()); // atom z position
		  atom_count++;
		  atom_list.push_back(atom_temp);          // push the atom
		  poscount=0;
		  break;
		default:
		  cout << "Something went terribly wrong" << endl;
		  return 0;
		  break;	  
		}
	    }
	}
      // Computing data
      if ( n == n_atoms*NB_ELEMENTS_LINE + BUFF_ELEMENTS )
	{
	  molecules=make_molecules(atom_list, cut_off_radius,box); // Making molecules from atom list
	  if( step == 0 ) // Bookeeping original molecules for diffusion 
	    {
	      molecules_legacy=molecules;
	    }
	  for( int i=0; i < molecules_legacy.size() ; i++)
	    {
	      for( int j=1; j < molecules_legacy[i].list.size(); j++)
		{
		  double d=distance_atoms(atom_list,molecules_legacy[i].list[0].index,molecules_legacy[i].list[j].index,box); // distance C-O for C and O in the same molecule at t=0
		  if ( d > cut_off_radius ) // if d > cut-off then add the count for diffused.
		    {
		      count_diff++;
		    }
		  diff_dist += d; // adding distance to the other for average
		  count_bond++;   // incrementing number of bonds
		}
	    }
	  // Diffusion
	  diff << step << " " << diff_dist/(double)count_bond << endl;
	  // Get different types of molecules
	  molecules1=getMoleculeWithN(molecules,2);
	  molecules2=getMoleculeWithN(molecules,3);
	  molecules3=getMoleculeWithN(molecules,4);
	  molecules4=getMoleculeWithN(molecules,5);
	  molecules5=getMoleculeWithN(molecules,6);
	  // Output coordinance values
	  coordC << step << " ";
	  coordC << molecules.size() << " ";
	  coordC << molecules1.size() << " ";
	  coordC<< molecules2.size() << " ";
	  coordC << molecules3.size() << " ";
	  coordC << molecules4.size() << " ";
	  coordC << molecules5.size() << endl;
	  // Distances count for different molecules
	  avg_dist << step << " ";
	  avg_dist << distance_avg_mols(molecules2,box) << " ";
	  avg_dist << distance_avg_mols(molecules3,box) << " ";
	  avg_dist << distance_avg_mols(molecules4,box) << " ";
	  avg_dist << distance_avg_mols(molecules,box)  << endl;
	  //====================================================================
	  firstNN_C << step << " ";
	  for( int i=0 ; i < 32 ; i++ )
	    {
	      firstNN_C << getNN_distance(atom_list,i,box) <<  " ";
	    }
	  firstNN_C << endl;
	  //===================================================================
	  firstNN_O << step << " " ;
	  for( int i=32 ; i < atom_list.size() ; i++ )
	    {
	      firstNN_O << getNN_distance(atom_list,i,box) <<  " ";
	    }
	  firstNN_O << endl;
	  //===================================================================
	  secondNN_C << step << " ";
	  for( int i=0 ; i < 32 ; i++ )
	    {
	      secondNN_C << getNN2_distance(atom_list,i,box) <<  " ";
	    }
	  secondNN_C << endl;
	  //==================================================================
	  secondNN_O << step << " " ;
	  //
	  for( int i=32 ; i < atom_list.size() ; i++ )
	    {
	      secondNN_O << getNN2_distance(atom_list,i,box) <<  " ";
	    }
	  secondNN_O << endl;
	  //==================================================================
	  // Iteration stuff
	  diff_dist  = 0; // Reinitialize diffusion distance
	  count_bond=0;  // number of bonds
	  count_diff=0;  // number of otters
	  atom_count=0; // Reinitialize count for atoms
	  atom_list.clear(); // Clear the list of atoms, for next step
	  step++;  // Incrementing step count
	  n=0;     // Reinitialize the element number for reading
	}
      ++it; // Reads next element
    }
  // Closing fluxes
  avg_dist.close();
  coordC.close();
  coordO.close();
  firstNN_C.close();
  firstNN_O.close();
  secondNN_C.close();
  secondNN_O.close();
  diff.close();
  // END 
  return 0;
}
