#include "lut.h"

//=======
// TESTS
//======================================================================
bool testType( TypeLUT lut , std::string specie )
{
  if ( lut.type_name == specie ) return true;
  else return false;
}
//---------------------------------------------------------------
bool testType( TypeLUT lut , int index )
{
  if ( lut.type_index == index ) return true;
  else return false;
}
//---------------------------------------------------------------
bool typeExists( std::vector<TypeLUT> list , std::string specie )
{
  for ( int i = 0 ; i < list.size() ; i++ )
    {
      if ( testType( list[i], specie ) ) return true;
    }
  return false;
}
//---------------------------------------------------------------
bool typeExists( std::vector<TypeLUT> list , int index )
{
  for ( int i = 0 ; i < list.size() ; i++ )
    {
      if ( testType( list[i], index ) ) return true;
    }
  return false;
}
//======================================================================

//==================
// ADD ATOM TO LUT
//======================================================================
void addAtom2LUT( std::vector<TypeLUT>  & list , std::string name , int index  )
{
  for ( int i=0 ; i < list.size() ; i++ )
    {
      if ( name == list[i].type_name )
	{
	  list[i].atom_index.push_back( index );
	  return;
	}
    }
  std::vector<int> atom_index;
  atom_index.push_back( index );
  TypeLUT lut = { name, list.size() , atom_index };
  list.push_back( lut );
  return;
}
//--------------------------------------------------------------
void addAtom2LUT( std::vector<TypeLUT> & list , Atom atom )
{
  for ( int i=0 ; i < list.size() ; i++ )
    {
      if ( atom.name == list[i].type_name )
	{
	  list[i].atom_index.push_back( atom.index );
	  return;
	}
    }
  std::vector<int> atom_index;
  atom_index.push_back(atom.index);
  TypeLUT lut = { atom.name, list.size() , atom_index };
  return ;
}
//======================================================================

//==========
// MAKE LUT
//======================================================================
TypeLUT makeLUT( std::string name , int index )
{
  std::vector<int> atom_index;
  return { name, index, atom_index };
}
TypeLUT makeLUT( std::string name , int index , std::vector<int> atom_index )
{
  return { name , index , atom_index };
}
//--------------------------------------------------------------
std::vector<TypeLUT> makeLUT( std::vector<Atom> atoms )
{
  std::vector<TypeLUT> list_lut;
  for ( int i=0 ; i < atoms.size() ; i++ )
    {
      addAtom2LUT( list_lut, atoms[i] );
    }
  return list_lut;
}
//======================================================================
