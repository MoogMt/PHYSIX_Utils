#ifndef CONTACT_MATRIX_H
#define CONTACT_MATRIX_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

#include "utils.h"
#include "atom.h"
#include "cell.h"

//================
// CONTACT MATRIX
//====================================================================
// Restricted Contact Matrix
struct Contact_Matrix
{
  std::vector<std::string> types;  // Types of the atoms
  std::vector<double> matrix;      // Contact matrix
};
//---------------------------------------------------------
// Full Contact Matrix
template <int nb_atoms>
struct ContactMatrix
{
  int types[nb_atoms];           // Atom types
  int matrix[nb_atoms*nb_atoms]; // Contact_Matrix
};
//---------------------------------------------------------
// Cut-Off Matrix
template <int n_type>
struct CutOffMatrix
{
  std::string types1[n_type]; std::string types2[n_type];
  double cut_offs_matrix[n_type];
};
//====================================================================

//===========
// FUNCTIONS
//====================================================================
// MAKE
//--------------------------------------------------------------------------------
// Computes the contact matrix from an atom list and a cell
Contact_Matrix makeContactMatrix ( std::vector<Atom> atom_list, Cell box );
template <int N> ContactMatrix<N> makeContactMatrix ( std::vector<Atom> atom_list, Cell box );
//--------------------------------------------------------------------------------
// ANGLES ET DISTANCE
//------------------------------------------------------------------------------------
double getAngle( Contact_Matrix contact_matrix, int atom_center_index , int atom_2_index, int atom_3_index );
double getDistance(Contact_Matrix contact_matrix, int atom_index_1, int atom_index_2 );
void writeAtomDistances( std::ofstream & file , std::vector<Atom> atom_list , std::vector<int> atom_index, Cell box);
//-------------------------------------------------------------------------------------
// CONTACT
//---------------------------------------------------------------------------------------
std::vector<double> getAtomContact(Contact_Matrix contact_matrix, int atom_index);
std::vector<double> getAtomContact( Contact_Matrix contact_matrix , int atom_index, std::string specie );
void writeAtomContact( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> atom_index );
//----------------------------------------------------------------------------------------
// COORDINANCE
//-----------------------------------------------------------------------------------------------
int getAtomNeighboursNb( Contact_Matrix contact_matrix, int atom_index, double cut_off_radius );
int getAtomNeighboursNb( Contact_Matrix contact_matrix, int atom_index, std::string specie, double cut_off_radius );
std::vector<int> getAtomsNeighboursNb( Contact_Matrix contact_matrix , std::vector<int> atom_index_list , double cut_off_radius );
std::vector<int> getTypeNeighboursNb(Contact_Matrix contact_matrix, std::string type, double cut_off_radius ) ;
double getTypeCoordinance( Contact_Matrix contact_matrix, std::string type, double cut_off_radius );
// Nearest Neighbours
//-------------------------------------------------------------------------------------------------
// -> Get
double getNNearest( Contact_Matrix contact_matrix , int n_nearest, int atom_index );
double getNNearest( Contact_Matrix contact_matrix , int n_nearest, int atom_index, std::string specie );
std::vector<double> getNNearest( Contact_Matrix contact_matrix , std::vector<int> n_nearest, int atom_index);
std::vector<double> getNNearest( Contact_Matrix contact_matrix , int nearest, std::vector<int> atom_indexes );
std::vector<double> getNNearest( Contact_Matrix contact_matrix , int nearest, std::vector<int> atom_indexes , std::string specie );
//-------------------------------------------------------------------------------
// TO BE REPLACED 
std::vector<std::vector<double> > getNNearest( Contact_Matrix contact_matrix , std::vector<int> nearest, std::vector<int> atom_indexes , std::string specie );
//-------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix contact_matrix , int n_nearest, std::string atom_type );
std::vector<double> getNNearest( Contact_Matrix contact_matrix , int n_nearest, std::vector<std::string> atom_types );
// -> Write
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , int n_nearest , std::vector<int> atom_indexes , int step );
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> nearest, int atom_index);
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> nearest, std::vector<int> atom_indexes , int step);
//====================================================================

#endif 
