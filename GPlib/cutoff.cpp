#include "cutoff.h"

//======
// READ
//====================================================================================
CutOffMatrix readCutOff( const std::string file , std::vector<TypeLUT> & lut_list )
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
      lut_list.push_back( makeLUT( names[i] , i ) );
    }
  //----------------------------------------------------------

  //-----------------
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
double getCutOff( const CutOffMatrix cut_off_matrix , const int i , const int j )
{
  return cut_off_matrix.matrix[ i*( cut_off_matrix.list_lut.size() ) + j ];
}
//============================================================================================
