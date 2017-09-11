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
};
//------------------------------

// molecules
//--------------------
struct molecule
{
  vector<atom> list;
};
//-------------------


// Print for debug
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
//--------------------------------------------
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
double pbc(double x1, double x2, double a)
{
  vector<double> x;
  x.push_back( pow( x1-x2   ,2) );
  x.push_back( pow( x1+a-x2 ,2) );
  x.push_back( pow( x1-a-x2 ,2) );
  return min(x);
}

// Distance
//------------------------------------------------------
double distance2(vector<atom> atoms, int i, int j, double a, double b, double c)
{
  double x=pbc(atoms[i].x,atoms[j].x,a);
  double y=pbc(atoms[i].y,atoms[j].y,b);
  double z=pbc(atoms[i].z,atoms[j].z,c);
  return sqrt(x+y+z);
}
//------------------------------------------------------


// Average
//---------------------------------------
double average(vector<double> data)
{
  double average=0;
  for( int i=0; i<data.size(); i++)
    {
      average=average+data[i];
    }
  return average=average/data.size();
}
//--------------------------------------

// Make molecules
//-------------------------------------------------------
vector<molecule> make_molecule( vector<atom> atom_list)
{
  // Variables
  vector<molecule> molecules;
  molecule molecule_buff;
  // Loop over atoms
  for ( int i=0 ; i < atom_list.size() ; i++ )
    {
      if ( (atom_list[i].name).compare("C") == 0 )
	{
	  (molecule_buff.list).push_back(atom_list[i]);
	  for ( int j=0 ; j < atom_list.size() ; j++)
	    {
	      if ( (atom_list[j].name).compare("O") == 0 && distance2(atom_list,i,j,9.0,9.0,9.0) < 1.75 && i!=j )
		{
		  (molecule_buff.list).push_back(atom_list[j]);
		}
	    }
	  molecules.push_back(molecule_buff);
	  (molecule_buff.list).clear();
	}
    }
  return molecules;
}
//-----------------------------------------------------


// GetMols
//--------------------------------------------------------
vector<molecule> getCOn(vector<molecule> mols, int n)
{
  vector<molecule> mols2;
  for( int i=0; i<mols.size(); i++)
    {
      if( (mols[i].list).size() == n)
	{
	  mols2.push_back(mols[i]);
	}
    }
  return mols2;
}
//--------------------------------------------------------

double distance_avg(molecule mol, int n)
{
  double distance=0;
  for ( int i=1; i < (mol.list).size() ; i ++ )
    {
      distance=distance + distance2(mol.list,0,2,9.0,9.0,9.0);
    }
  return distance/(double)(n);
}

double distance_avg(vector<molecule> mols, int n)
{
  double distance=0;
  if ( mols.size() == 0 )
    {
      return 0;
    }
  for( int i=0; i < mols.size() ; i++)
    {      distance = distance + distance_avg(mols[i], n);
    }
  return distance/((double)(mols.size()));
}

double distance_avg(vector<molecule> mols)
{
  double distance=0;
  int count=0;
  for (int i=0; i<mols.size(); i++)
    {
      distance = distance + distance_avg(mols[i],(mols[i].list).size()-1);
      count++;
    }
  return distance/(double)count;
}

//=====================================================================
int main(void)
{
  // Tab
  //-----------------------
  vector<atom> atom_list;
  atom atom_buff;
  vector<molecule> mols;
  vector<molecule> mols2;
  vector<molecule> mols3;
  vector<molecule> mols4;
  
  //-----------------------
  
  // Input
  //------------------------------------
  ifstream fichier("traj.xyz");
  // Input Flux Iterator
  istream_iterator<string> it(fichier);
  istream_iterator<string> end;
  //------------------------------------
  
  // Count
  //--------------------
  int n=0;
  int step=0;
  int poscount=0;
  //-------------------

  // Data
  vector<double> angles;
  vector<double> distances;
  vector<double> coordC;
  vector<double> coordO;

  // Reading the file 
  while(it != end )   
    {
      //------------------------
      // Item number iteration
      //------------------------
      ++n;
      //------------------------
      if ( n > 4 )
	{
	  if( string(*it).compare("C") == 0 ||  string(*it).compare("O") == 0 )
	    {
	      atom_buff.name=string(*it);
	    }
	  else
	    {
	      poscount++;
	      switch(poscount)
		{
		case 1:
		  atom_buff.x = atof(string(*it).c_str());
		  break;
		case 2:
		  atom_buff.y = atof(string(*it).c_str());
		  break;
		case 3:
		  atom_buff.z = atof(string(*it).c_str());
		  atom_list.push_back(atom_buff);
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
      if ( n == 32*3*4+3 )
	{
	  // Making molecules
	  //------------------------------
	  mols=make_molecule(atom_list);
	  //------------------------------
	  // Crunching data
	  //---------------------------------------------------
	  mols2=getCOn(mols,3);
	  mols3=getCOn(mols,4);
	  mols4=getCOn(mols,5);
	  //cout << step << " " << mols.size() << " " << mols2.size() << " " << mols3.size() << " " << mols4.size() << endl;
	  //cout << step << " " <<  distance_avg(mols2,2) << " " << distance_avg(mols3,3) << " " << distance_avg(mols4,4) << endl;
	  if ( distance_avg(mols4,4) != 0 )
	    {
	      cout << step << " " <<  distance_avg(mols2,2) << " " << distance_avg(mols3,3) << " " <<  distance_avg(mols) << " " << distance_avg(mols4,4)  << endl;
	    }
	  else
	    {
	      cout << step << " " <<  distance_avg(mols2,2) << " " << distance_avg(mols3,3) << " " <<  distance_avg(mols) << endl; 
	    }
	  //---------------------------------------------------
	  // Clear
	  //------------------
	  atom_list.clear();
	  //------------------
	  // Boost
	  //-------------------
	  step++;  // Next step
	  n=0;     // reinitialize atoms
	  //-------------------
	}
      //-------------------------------------------------------
      //------------------
      // Go to next item
      //------------------
      ++it;
      //------------------
    }
  //---------------
  // Data averages
  //-------------------------

  //-------------------------
  // Printing data
  //---------------
  //---------------
  return 0;
}
