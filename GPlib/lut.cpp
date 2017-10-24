#include "lut.h"

//=======
// TESTS
//======================================================================
bool testType ( const AllTypeLut lut_all, const int atom_index ,  const std::string specie )
{
  if ( lut_all.type_name[ atom_index ] == specie ) ) return true;
  else return false;
}
//---------------------------------------------------------------
bool testType ( const AllTypeLut lut_all, const int atom_index ,  const int type_index )
{
  if ( lut_all.type_index[ atom_index ] == type_index ) ) return true;
  else return false;
}
//---------------------------------------------------------------
bool testType( const TypeLUT lut , const std::string specie )
{
  if ( lut.type_name == specie ) return true;
  else return false;
}
//---------------------------------------------------------------
bool testType( const TypeLUT lut , const int index )
{
  if ( lut.type_index == index ) return true;
  else return false;
}
//---------------------------------------------------------------
bool typeExists( const std::vector<TypeLUT> list , const std::string specie )
{
  for ( int i = 0 ; i < list.size() ; i++ )
    {
      if ( testType( list[i], specie ) ) return true;
    }
  return false;
}
//---------------------------------------------------------------
bool typeExists( const std::vector<TypeLUT> list , const int index )
{
  for ( int i = 0 ; i < list.size() ; i++ )
    {
      if ( testType( list[i], index ) ) return true;
    }
  return false;
}
//---------------------------------------------------------------
bool typeExists( const AllTypeLUT lut , const std::string specie )
{
  if ( typeExists( lut.types, specie ) ) return true;
  else return false;
}
//---------------------------------------------------------------
bool typeExists( const AllTypeLUT lut , const int type_index )
{
  if ( typeExists( lut.types, type_index ) ) return true;
  else return false;
}
//======================================================================

//==================
// ADD ATOM TO LUT
//======================================================================
void addAtom2LUT( std::vector<TypeLUT>  & list , const std::string name , const int index  )
{
  for ( int i=0 ; i < list.size() ; i++ )
    {
      if ( name == list[i].type_name )
	{
	  list[i].atom_index.push_back( index );
	  return;
	}
    }
  list.push_back( { name, list.size() , initVector( index ) } );
  return;
}
//--------------------------------------------------------------
void addAtom2LUT( std::vector<TypeLUT> & list , const Atom atom )
{
  for ( int i=0 ; i < list.size() ; i++ )
    {
      if ( atom.name == list[i].type_name )
	{
	  list[i].atom_index.push_back( atom.index );
	  return;
	}
    }
  list.push_back( { atom.name , list.size() , initVector( atom.index ) } );
  return ;
}
//--------------------------------------------------------------
void addAtom2LUT( AllTypeLUT & list , const std::string name , const int index  )
{
  list.type_name.push_back( name );
  list.type_index.push_back( index );
  addAtom2LUT( list.types , name , index );
}
//--------------------------------------------------------------
void addAtom2LUT( AllTypeLUT & list , const Atom atom )
{
  list.type_name.push_back( atom.name );
  list.type_index.push_back( atom.index );
  addAtom2LUT( list.types , atom.name , atom.index );
} 
//======================================================================

//==========
// MAKE LUT
//======================================================================
TypeLUT makeLUT( const std::string name , const int index )
{
  std::vector<int> atom_index;
  return { name, index, atom_index };
}
//--------------------------------------------------------------
TypeLUT makeLUT( const std::string name , const int index , const std::vector<int> atom_index )
{
  return { name , index , atom_index };
}
//--------------------------------------------------------------
std::vector<TypeLUT> makeLUT( const std::vector<Atom> atoms )
{
  std::vector<TypeLUT> list_lut;
  for ( int i=0 ; i < atoms.size() ; i++ )
    {
      addAtom2LUT( list_lut, atoms[i] );
    }
  return list_lut;
}
//======================================================================