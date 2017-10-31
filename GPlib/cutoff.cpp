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
	  if ( !( it_string >> com.matrix[i*nb_type+j] ) ) return false;
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
  //--------------------------
  // Stream related variables
  //-------------------------------------------
  std::ifstream input( file.c_str() );

  if ( readCutOff( input , com , lut_list ) ) return true;
  else
    {
      std::cout << "Problem reading file " << file << std::endl;
      return false;
    }
}
//---------------------------------------------------------------------------------------
CutOffMatrix readCutOff( const std::string file , AllTypeLUT & lut_list )
{
  //-----------
  // Variables
  //-----------------
  CutOffMatrix com;
  //-----------------
  
  //--------------------------
  // Stream related variables
  //-------------------------------------------
  std::ifstream input( file.c_str() );
  std::istream_iterator<std::string> read(input);
  std::istream_iterator<std::string> end;
  //-------------------------------------------

  ///----------
  // Variables
  //--------------------------------------------------------------
  int count = 0; // count
  int nb_types = atoi( std::string( *read ).c_str() ); ++read; // Reads number of types
  //--------------------------------------------------------------

  //------------------------
  // Reading names of types
  //----------------------------------------------------------
  std::vector<std::string> names;
  while( read != end && count < nb_types )
    {
      names.push_back( std::string( *read ) );
      ++read;
      count++;
    }
  //----------------------------------------------------------
  
  //----------------
  // BUILDING LUT
  //----------------------------------------------------------
  for ( int i=0 ; i < nb_types ; i++ )
    {
      lut_list.types.push_back( makeLUT( names[i] , i ) );
    }
  //----------------------------------------------------------

  //------------------
  // READING MATRIX
  //----------------------------------------------------------
  int i=0;
  while ( read != end && i < nb_types*nb_types )
    {
      com.matrix.push_back( atof( std::string(*read).c_str() ) ); ++read;
      i++;
    }
  //----------------------------------------------------------
  
  return com;
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
