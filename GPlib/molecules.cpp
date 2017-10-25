#include "molecules.h"

//================
// Make molecules
//=================================================================================
MoleculeBasic startMolecule( const std::string name_atom , const int index_atom )
{
  std::vector<std::string> names;
  std::vector<int> index;
  names.push_back( name_atom );
  index.push_back( index_atom );
  return { names, index };
}
//----------------------------------------------------------------------------
std::vector<MoleculeBasic> makeMolecules( const ContactMatrix & cm )
{
  //-----------------
  // Initialization
  //------------------------------------------------------------------
  std::vector<MoleculeBasic> mol_list;
  int nb_atoms = (int)( sqrt( cm.matrix.size() ) );
  // Initiate all connection vectors
  int used[ nb_atoms ]; zeros( used , nb_atoms ); // check if atom is already connected
  int try1[ nb_atoms ]; zeros( try1 , nb_atoms ); // global connection vector
  int try2[ nb_atoms ]; zeros( try2 , nb_atoms ); // local atomic connection vector
  //------------------------------------------------------------------

  //-----------------------------------------------------------------------
  int i = 0; // dummy index
  // Loop over all non already used atoms
  while ( i < nb_atoms && used[i] == 0 )
    {
      // Reinitiate matrix
      zeros( try1 , nb_atoms ); try1[i] = 1;
      if ( cm.lut_list.type_name.size() <= i ) exit(0);
      MoleculeBasic molecule = startMolecule( cm.lut_list.type_name[i] , cm.lut_list.type_index[i] );
      // Starting molecule
      do
	{
	  // Reinitiate matrix
	  zeros( try2 , nb_atoms ) ;
	  // Second loop over all atoms
	  for ( int k=0 ; k < nb_atoms ; k++ )
	    {
	      if ( try1[k] == 0 )
		{
		  continue;
		}
	      used[k] = 1;
	      // Loop over all atoms
	      for ( int h=0 ; h < nb_atoms ; h++ )
		{
		  if ( used[h] == 0 && connected(cm,k,h) )
		    {
		      try2[h] = 1;
		      molecule.names.push_back( cm.lut_list.type_name[h] );
		      molecule.atom_index.push_back( cm.lut_list.type_index[h] );
		    }
		}
	    }
	  // Propagates the neighbor
	  copy( try2 , try1 , nb_atoms);
	} while ( sum( try1, nb_atoms ) == 0 );
      // Add molecule to list
      mol_list.push_back( molecule );
      // Going to the next atom
      i++;
    }
  //-----------------------------------------------------------------------
  return mol_list;
}
//=================================================================================

//================
// PRINT MOLECULE
//=================================================================================
void printMolecule( const MoleculeBasic molecule )
{
  for ( int i=0 ; i < molecule.names.size() ; i++ )
    {
      std::cout << molecule.names[i] << " " << molecule.atom_index[i] << " ";
    }
  std::cout << std::endl;
  return;
}
//-----------------------------------------------------------------------------
void printMoleculeSize( const MoleculeBasic molecule, const bool toline )
{
  std::cout << molecule.names.size() << " ";
  if ( toline ) std::cout << std::endl;
  return;
}
//-----------------------------------------------------------------------------
void printMolecules( const std::vector<MoleculeBasic> molecules )
{
  std::cout << "------------------------" << std::endl;
  for ( int i=0 ; i < molecules.size() ; i++ )
    {
      printMolecule( molecules[i] );
    }
  std::cout << "------------------------" << std::endl;
  return ;
}
//-----------------------------------------------------------------------------
void printMoleculesSize( const std::vector<MoleculeBasic> molecules )
{
  for ( int i=0 ; i < molecules.size() ; i++ )
    {
      printMoleculeSize( molecules[i], false );
    }
  std::cout << std::endl;
  return;
}
//=================================================================================
