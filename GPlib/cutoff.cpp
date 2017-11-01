#include "cutoff.h"

//======
// READ
//====================================================================================
bool readCutOff( std::ifstream & input , CutOffMatrix & com , std::vector<TypeLUT> & lut_list )
{
  //---------------
  // String methods
  //----------------------------------------------
  std::string line;
  //----------------------------------------------

  //----------------------------
  // Getting the number of types
  //-----------------------------------------------
  int nb_type;
  if ( std::getline( input , line ) )
    {
      std::istringstream it_string(line);
      if ( ! ( it_string >> nb_type ) ) return false;
    }
  //-----------------------------------------------

  //---------------------
  // Reading types names
  //-----------------------------------------------
  int i=0;
  std::vector<std::string> names; names.assign( nb_type, "" );
  while( i < nb_type && getline( input, line ) )
    {
      std::istringstream it_string(line);
      std::string name;
      if ( ! ( it_string >> name) ) return false;
      names[i] = name ;
      i++;
    }
  //-----------------------------------------------
  
  //------------------------
  // Reading names of types
  //----------------------------------------------------------
  //  std::vector<std::string> names;
  /*while( read != end && count < nb_types )
    {
      names.push_back( std::string( *read ) );
      ++read;
      count++;
      }*/
  //----------------------------------------------------------

  // Checking names
  if ( !checkNames(names) ) return false;

  // BUILDING LUT
  makeLUT( lut_list , names);

  // Reading matrix
  i=0;
  com.matrix.assign( nb_type*nb_type, 0.0);
  while( getline( input, line) && i < nb_type )
    {
      int j=i;
      std::istringstream it_string(line);
      while( j < nb_type )
	{
	  if ( !( it_string >> com.matrix[ i*nb_type+j ] ) ) return false;
	  com.matrix[ j*nb_type + i ];
	  j++;
	}
    }

  /*while ( read != end && i < nb_types*nb_types )
    {
      com.matrix.push_back( atof( std::string(*read).c_str() ) ); ++read;
      i++;
      }*/
  
  return true;
}
//--------------------------------------------------------------------------------------
bool readCutOff( const std::string file , CutOffMatrix & com , std::vector<TypeLUT> & lut_list )
{
  std::ifstream input( file.c_str() );

  if ( readCutOff( input , com , lut_list ) ) return true;
  else
    {
      std::cout << "Problem reading file " << file << std::endl;
      return false;
    }
}
//---------------------------------------------------------------------------------------
bool readCutOff(  std::ifstream & input , CutOffMatrix & com, AllTypeLUT & lut_list )
{
    //---------------
  // String methods
  //----------------------------------------------
  std::string line;
  //----------------------------------------------

  //----------------------------
  // Getting the number of types
  //-----------------------------------------------
  int nb_type;
  if ( std::getline( input , line ) )
    {
      std::istringstream it_string(line);
      if ( ! ( it_string >> nb_type ) ) return false;
    }
  //-----------------------------------------------

  //---------------------
  // Reading types names
  //-----------------------------------------------
  int i=0;
  std::vector<std::string> names; names.assign( nb_type, "" );
  while( i < nb_type && getline( input, line ) )
    {
      std::istringstream it_string(line);
      std::string name;
      if ( ! ( it_string >> name) ) return false;
      names[i] = name ;
      i++;
    }
  //-----------------------------------------------

  //----------
  // Build LUT
  //----------------------------------------------
  // Checking Names
  if ( !checkNames(names) ) return false;
  // Building LUT
  makeLUT( lut_list , names);
  //----------------------------------------------

  //----------------
  // Reading matrix
  //----------------------------------------------
  i=0;
  // Assign Matrix
  com.matrix.assign( nb_type*nb_type, 0.0);
  // Reading matrix
  while( getline( input, line) && i < nb_type )
    {
      int j=i;
      std::istringstream it_string(line);
      while( j < nb_type )
	{
	  if ( !( it_string >> com.matrix[i*nb_type+j] ) ) return false;
	  com.matrix[ j*nb_type + i ];
	  j++;
	}
    }
  //----------------------------------------------

  return true;
}
//-----------------------------------------------------------------------------------
bool readCutOff( const std::string file , CutOffMatrix & com, AllTypeLUT & lut_list )
{
  // File Handler
  std::ifstream input( file.c_str() );
  // Reading Cut_Off File
  if ( readCutOff( input, com, lut_list ) ) return true;
  // If problem...
  else
    {
      std::cout << "Problem reading file: " << file << std::endl;
      return false;
    }
}
//====================================================================================

//=====
// Get
//====================================================================================
double getCutOff( const CutOffMatrix cut_off_matrix , const int type_index1 , const int type_index2 )
{
  int nb_type = sqrt( cut_off_matrix.matrix.size() );
  return cut_off_matrix.matrix[ type_index1*nb_type + type_index2 ];
}
//============================================================================================
