#include "utils.h"

//======
// MIN
//===========================================================
int min( int int1, int int2 )
{
  if (int1 < int2) return int1;
  else return int2;
}
//-------------------------------------------
double min ( double real1, double real2 )
{
  if ( real1 < real2 ) return real1;
  else return real2;
}
//-------------------------------------------
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
//===========================================================

//=====
// MAX
//===========================================================
int max( int int1, int int2 )
{
  if (int1 > int2) return int1;
  else return int2;
}
//-------------------------------------------
double max ( double real1, double real2 )
{
  if ( real1 > real2 ) return real1;
  else return real2;
}
//-------------------------------------------
double max(std::vector<double> vector)
{
  double max=vector[0];
  for( int i=0 ; i < vector.size() ; i++ )
    {
      if( max < vector[i] )
	{
	  max=vector[i];
	}
    }
  return max;
}
//==========================================================

//=========
// AVERAGE
//==========================================================
double average(std::vector<int> data)
{
  double average=0;
  for( int i=0; i<data.size(); i++)   
    {
      average=average+(double)(data[i]);
    }
  return average=average/(double)(data.size()); // Returning the average
}
//----------------------------------------------------------
double average(std::vector<double> data)
{
  double average=0;                   // Average
  for( int i=0; i<data.size(); i++)   
    {
      average=average+data[i];
    }
  return average=average/data.size(); // Returning the average
}
//==========================================================

//======
// SUMS
//==========================================================
int sumFromO(int integer)
{
  int sum=0;
  for ( int i=0 ; i <= integer ; i++ )
    {
      sum += i;
    }
  return sum;
}
//-------------------------------------------------------
int sumBtw(int int1, int int2)
{
  int sum=0;
  for ( int i=min(int1,int2) ; i < max(int1,int2) ; i++ )
    {
      sum += i;
    }
  return sum;
}
//-------------------------------------------------------
int computeSep(int atom_index, int nb_atoms)
{
  return atom_index*nb_atoms-sumFromO(atom_index);
}
//==========================================================

//========
// VECTORS
//==============================================================================
void switchV( std::vector<int> & vector , int index1 , int index2 )
{
  double stock = vector[index1];
  vector[ index1 ] = vector[ index2 ];
  vector[ index2 ] = stock; 
  return;
}
//-------------------------------------------------------
void switchV( std::vector<double> & vector , int index1 , int index2 )
{
  double stock = vector[index1];
  vector[ index1 ] = vector[ index2 ];
  vector[ index2 ] = stock; 
  return;
}
//-------------------------------------------------------
std::vector<int> initVector( int value )
{
  std::vector<int> vector;
  vector.push_back( value );
  return vector;
}
//-------------------------------------------------------
std::vector<double> initVector( double value )
{
  std::vector<double> vector;
  vector.push_back( value );
  return vector;
}
//-------------------------------------------------------
std::vector<std::string> initVector( std::string value )
{
  std::vector<std::string> vector;
  vector.push_back( value );
  return vector;
}
//-------------------------------------------------------
std::vector<int> makeVec( int init, int final )
{
  std::vector<int> vector;
  for ( int i = init ; i <= final ; i++ )
    {
      vector.push_back(i);
    }
  return vector;
}
//-----------------------------------------------------------------------
std::vector<double> sortVector(std::vector<double> to_sort, bool increasing)
{
  if ( increasing )
    {
      return sortVectorIncreasing(to_sort);
    }
  else
    {
      return sortVectorDecreasing(to_sort);
    }
}
//-----------------------------------------------------------------------
std::vector<double> sortVectorIncreasing(std::vector<double> to_sort)
{
  for ( int i=0 ; i < to_sort.size()-1 ; i++ )
    {
      for ( int j=i ; j < to_sort.size() ; j++ )
	{
	  if ( to_sort[i] > to_sort[j] )
	    {
	      double stock = to_sort[i];
	      to_sort[i] = to_sort[j];
	      to_sort[j] = stock;
	    }
	}
    }
  return to_sort;
}
//-----------------------------------------------------------------------
std::vector<double> sortVectorDecreasing( std::vector<double> to_sort )
{
  for ( int i=0 ; i < to_sort.size()-1 ; i++ )
    {
      for ( int j=i ; j < to_sort.size() ; j++ )
	{
	  if ( to_sort[i] < to_sort[j] )
	    {
	      double stock = to_sort[i];
	      to_sort[i] = to_sort[j];
	      to_sort[j] = stock;
	    }
	}
    }
  return to_sort;
}
//-----------------------------------------------------------------------
void appendVector( std::vector<double> & vec_to , std::vector<double> vec_from )
{
  for ( int i=0 ; i < vec_from.size() ; i++ )
    {
      vec_to.push_back( vec_from[i] );
    }
  return;
}
//-----------------------------------------------------------------------
std::vector<double> difference( std::vector<double> vector1 , std::vector<double> vector2 )
{
  if ( vector1.size() == vector2.size() )
    {
      for ( int i=0 ; i < vector1.size() ; i++ )
	{
	  vector1[i] -= vector2[i];
	}
      return vector1;
    }
  else vector1.clear();
}
//-----------------------------------------------------------------------
double average( std::vector<double> & vector )
{
  double avg = 0;
  for ( int i=0 ; vector.size() ; i++ )
    {
      avg += vector[i];
    }
  return avg/(double)(vector.size());
}
//-----------------------------------------------------------------------
std::vector<double> distance( std::vector<double> vector1 , double scalar )
{
  for ( int i=0 ; i < vector1.size() ; i++ )
    {
      vector1[i] -= scalar; 
    }
  return vector1;
}
//-----------------------------------------------------------------------
std::vector<double> square( std::vector<double> vector1 )
{
  for ( int i=0 ; i < vector1.size() ; i++ )
    {
      vector1[i] *= vector1[i];
    }
  return vector1;
}
//-----------------------------------------------------------------------
std::vector<double> squaroot( std::vector<double> vector1 )
{
  for ( int i=0; i < vector1.size() ; i++ )
    {
      vector1[i] = sqrt( vector1[i] );
    }
  return vector1;
}
//-----------------------------------------------------------------------
std::vector<double> addVector( std::vector<double> vector1 , std::vector<double> vector2 )
{
  if( vector1.size() == vector2.size() )
    {
      for( int i=0 ; i < vector1.size() ; i++ )
	{
	  vector1[i] += vector2[i];
	}
      return vector1;
    }
  else vector1.clear();
}
//-----------------------------------------------------------------------
double scalarProduct( std::vector<double> vector1 , std::vector<double> vector2 )
{
  double scalar = 0;
  if ( vector1.size() == vector2.size() )
    {
      for ( int i=0 ; i < vector1.size() ; i++ )
	{
	  scalar += vector1[i]*vector2[i];
	}
      return scalar;
    }
  else return 0;
}
//-----------------------------------------------------------------------
std::vector<double> crossProduct( std::vector<double> vector1 , std::vector<double> vector2 )
{
  std::vector<double> vector;
  vector.push_back( vector1[1]*vector2[2] - vector2[1]*vector1[2] );
  vector.push_back( vector1[2]*vector2[0] - vector2[2]*vector1[0] );
  vector.push_back( vector1[0]*vector2[1] - vector2[0]*vector1[1] );
  return vector;
}
//-----------------------------------------------------------------------
double norm( std::vector<double> vector )
{
  return sqrt(abs( scalarProduct( vector , vector ) ));
}
//-----------------------------------------------------------------------
double getDistanceFromPlan( std::vector<double> vector1, std::vector<double> vector2 , std::vector<double> point_outside , std::vector<double> point_plan )
{
  std::vector<double> normal = crossProduct( vector1, vector2 );
  std::vector<double> vector = difference( point_plan , point_outside  );
  return abs(scalarProduct( normal , vector ) ) / norm( normal );
}
//=================================================================

//=============
// CONVERSION
//=================================================================
double it2real (  std::istream_iterator<std::string> iterator )
{
  return atof(std::string(*iterator).c_str());
}
//=================================================================

//=======
// ARRAY
//======================================================
void zeros( int* vec , int nb_atoms )
{
  for ( int i=0 ; i < nb_atoms ; i++ )
    {
      vec[i] = 0;
    }
  return;
}
//------------------------------------------------
int sum( int* vector , int nb_atoms )
{
  double sum=0;
  for ( int i=0 ; i < nb_atoms ; i++ )
    {
      sum += vector[i];
    }
  return sum;
}
//------------------------------------------------
void copy( int* vec_from , int* vec_to , int size )
{
  for ( int i=0 ; i < size ; i++ )
    {
      vec_to[i] = vec_from[i];
    }
  return;
}
//======================================================
