#include "utils.h"

// Min
// Returns the minimum of value of vector containing doubles
//-----------------------------------------------------------
double min(std::vector<double> vector)
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

// Average
// Returns the average distance between two atoms
//----------------------------------------------------------
double average(std::vector<int> data)
{
  double average=0;
  for( int i=0; i<data.size(); i++)   
    {
      average=average+(double)(data[i]);
    }
  return average=average/(double)(data.size()); // Returning the average
}

double average(std::vector<double> data)
{
  double average=0;                   // Average
  for( int i=0; i<data.size(); i++)   
    {
      average=average+data[i];
    }
  return average=average/data.size(); // Returning the average
}
//--------------------------------------

// SUMS
int sumFromO(int integer)
{
  int sum=0;
  for ( int i=0 ; i <= integer ; i++ )
    {
      sum += i;
    }
  return sum;
}
int sumBtw(int int1, int int2)
{
  int sum=0;
  for ( int i=min(int1,int2) ; i < max(int1,int2) ; i++ )
    {
      sum += i;
    }
  return sum;
}


int computeSep(int atom_index, int nb_atoms)
{
  return atom_index*nb_atoms-sumFromO(atom_index);
}

int max( int int1, int int2 )
{
  if (int1 > int2) return int1;
  else return int2;
}

int min( int int1, int int2 )
{
  if (int1 < int2) return int1;
  else return int2;
}

std::vector<int> makeVec(int init, int final)
{
  std::vector<int> vector;
  for ( int i = init ; i < final ; i++ )
    {
      vector.push_back(i);
    }
  return vector;
}
