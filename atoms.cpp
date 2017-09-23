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
