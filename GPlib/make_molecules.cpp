//---------------
// External Libs
//-------------------
#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>
//---------------------

//----------------
// Internal Files
//-------------------------
#include "utils.h"
#include "atom.h"
#include "cell.h"
#include "contact_matrix.h"
#include "xyz.h"
#include "histogram.h"
#include "pdb.h"
#include "sim.h"
#include "lut.h"
#include "molecules.h"
//-------------------------

//================
// MAIN PROGRAM
//=====================================================================
int main( void )
{
  //--------
  // Input
  //---------------------------------
  std::ifstream input("TRAJEC.xyz");
  //---------------------------------

  //--------
  // Output
  //--------------------------------------------------------------------------------
  std::ofstream graph_molecules ( "graph_molecules.dat" ,  std::ios::out );
  std::ofstream hist_molecules ( "hist_molecules.dat" ,  std::ios::out );
  std::ofstream co2_stock_in_out ( "co2_stock_in_out.dat" ,  std::ios::out );
  std::ofstream co3_stock_in_out ( "co3_stock_in_out.dat" ,  std::ios::out );
  std::ofstream co4_stock_in_out ( "co4_stock_in_out.dat" ,  std::ios::out );
  std::ofstream co2_stock_in_time ( "co2_stock_in_time.dat" ,  std::ios::out );
  std::ofstream co3_stock_in_time ( "co3_stock_in_time.dat" ,  std::ios::out );
  std::ofstream co4_stock_in_time ( "co4_stock_in_time.dat" ,  std::ios::out );
  std::ofstream co2_stock_alone_out ( "co2_stock_alone_out.dat" ,  std::ios::out );
  std::ofstream co3_stock_alone_out ( "co3_stock_alone_out.dat" ,  std::ios::out );
  std::ofstream co4_stock_alone_out ( "co4_stock_alone_out.dat" ,  std::ios::out );
  std::ofstream co2_stock_alone_time ( "co2_stock_alone_time.dat" ,  std::ios::out );
  std::ofstream co3_stock_alone_time ( "co3_stock_alone_time.dat" ,  std::ios::out );
  std::ofstream co4_stock_alone_time ( "co4_stock_alone_time.dat" ,  std::ios::out );
  std::ofstream c2_ratio_out ( "c2_ratio_out.dat" ,  std::ios::out );
  std::ofstream c3_ratio_out ( "c3_ratio_out.dat" ,  std::ios::out );
  std::ofstream c4_ratio_out ( "c4_ratio_out.dat" ,  std::ios::out );
  std::ofstream distance_from_plan( "distance_from_plan.dat" , std::ios::out );
  std::ofstream distance_from_plan_in( "distance_from_plan_in.dat" , std::ios::out );
  std::ofstream distance_from_plan_alone( "distance_from_plan_alone.dat" , std::ios::out );
  std::ofstream distance_from_plan_hist_alone( "distance_from_plan_hist_alone.dat" , std::ios::out );
  std::ofstream distance_from_plan_hist_in( "distance_from_plan_hist_in.dat" , std::ios::out );
  //--------------------------------------------------------------------------------
  
  //----------------------
  // Physical parameters
  //--------------------------------------
  int step       = 1;  // Step counter
  int start_step = 2000; // Start step
  int end_step   = 1000000; // End step
  int comp_step  = 1; // Frequency of computation
  //--------------------------------------

  //---------------
  // Initializers
  //--------------------------------------------------
  AtomList  atom_list;  // Atoms in cell
  AllTypeLUT lut_list; // LUT for types
  ContactMatrix cm_connection;    // Contact Matrix
  ContactMatrix cm_distance;    // Contact Matrix
  //--------------------------------------------------

  //------
  // Data
  //-----------------------------------------
  std::vector<double> co2_data;
  std::vector<double> sizes;
  //-----------------------------------------

  //-----------
  // Histogram
  //-----------------------------------------
  // Sizes
  double hist_start = 0.5, hist_end = 96.5;
  double hist_start_ratio = 0. , hist_end_ratio = 1;
  int nb_box_ratio = 100;
  int nb_box = 96;
  std::vector<Bin> hist;
  //-----------------------------------------
  // Co2 %
  double hist_co2_start = 0.5, hist_co2_end = 32.5;
  int nb_co2_box = 32;
  std::vector<Bin> hist_co2;
  //-----------------------------------------
  int nb_box2 = 1000;
  //-----------------------------------------


  // CO2_fractions
  //-----------------------------------------
  int co2_in = 0, co2_alone = 0; std::vector<double> co2_stock_in; std::vector<double> co2_stock_alone;
  int co3_in = 0, co3_alone = 0; std::vector<double> co3_stock_in; std::vector<double> co3_stock_alone;
  int co4_in = 0, co4_alone = 0; std::vector<double> co4_stock_in; std::vector<double> co4_stock_alone;
  int o2_in = 0, o2_alone = 0; std::vector<double> o2_stock_in; std::vector<double> o2_stock_alone;
  int o3_in = 0, o3_alone = 0; std::vector<double> o3_stock_in; std::vector<double> o3_stock_alone;
  int c2_in = 0, c2_alone = 0; std::vector<double> c2_ratio ;
  int c3_in = 0, c3_alone = 0; std::vector<double> c3_ratio ;
  int c4_in = 0, c4_alone = 0; std::vector<double> c4_ratio ;
  std::vector<double> o2_ratio ;
  std::vector<double> o3_ratio ;
  std::vector<double> DistPlanC3;
  std::vector<double> distPlanC3_in ;
  std::vector<double> distPlanC3_alone ;
  //-----------------------------------------
  
  //--------------------
  // Reading Cell File
  //-------------------------------------------------------------------
  Cell cell;
  if ( ! readParamCell( "cell.param" , cell ) )
    {
      return 1;
    }
  //-------------------------------------------------------------------

  //-----------------
  // Reading Cut-Off
  //-------------------------------------------------------------------
  CutOffMatrix cut_off;
  if ( ! readCutOff( "cut_off.dat" , cut_off , lut_list ) )
    {
      return 1;
    }
  //-------------------------------------------------------------------

  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step % comp_step == 0 && step > start_step && step < end_step)
	{
       	  // Makes the contact matrix
	  makeContactMatrix( cm_connection , cm_distance , atom_list, cell , cut_off , lut_list );
	  // Making molecules
	  std::vector<Molecule> molecules = makeMolecules( cm_connection );

	  //----------------------------------------------------
	  // Printing bonds to file
	  //----------------------------------------------------
	  for( int i=0 ; i < molecules.size() ; i++ )
	    {
	      for ( int j=0 ; j < molecules[i].bonds.size(); j++ )
		{
		  int index_atom1 = molecules[i].bonds[j].atom1_index;
		  int index_atom2 = molecules[i].bonds[j].atom2_index;
		  graph_molecules << index_atom1 << " " << index_atom2 << " " << cm_connection.lut_list.type_index[ index_atom1 ] << " " << cm_connection.lut_list.type_index[ index_atom2 ] << std::endl;
		}
	      graph_molecules << "------------" << std::endl;
	    }
	  graph_molecules << "================" << std::endl;
	  //----------------------------------------------------

	  // Reinitialialize
	  co2_in = 0, co2_alone = 0;
	  co3_in = 0, co3_alone = 0;
	  co4_in = 0, co4_alone = 0;
	  
	  //--------------------
	  // Calculating angles
	  //-----------------------------------------------------------------------------
	  for ( int i=0 ; i < molecules.size() ; i++  )
	    {
	      for ( int j=0 ; j < molecules[i].atom_index.size() ; j++ )
		{
		  std::vector<double> angles = getAngleAtom( cm_distance , molecules[i] , molecules[i].atom_index[j] );
		  if ( angles.size() == 0 ) continue;
		  if ( molecules[i].names[j] == "C" )
		    {
		      if ( angles.size() == 1 )
			{
			  if ( molecules[j].names.size() == 3 )
			    {
			      co2_alone++;
			    }
			  else co2_in++;
			}
		      else if ( angles.size() == 3 )
			{
			  if ( molecules[j].names.size() == 4 )
			    {
			      co3_alone++;
			      std::vector<int> per_atoms_index = getBonded( molecules[i] , molecules[i].atom_index[j] );
			      int index_center = molecules[i].atom_index[j];
			      std::vector<double> position_center = getPosition( atom_list , index_center );
			      std::vector<double> position_atom1  = getPosition( atom_list , per_atoms_index[0] );
			      std::vector<double> position_atom2  = getPosition( atom_list , per_atoms_index[1] );
			      std::vector<double> position_atom3  = getPosition( atom_list , per_atoms_index[2] );
			      std::vector<double> vector1 = difference( position_atom1 , position_atom2 );
			      std::vector<double> vector2 = difference( position_atom1 , position_atom3 );
			      double dist = getDistanceFromPlan( vector1 , vector2 , position_center, position_atom1 );
			      distance_from_plan << step << " " << dist << std::endl;
			      DistPlanC3.push_back( dist );
			    }
			  else
			    {
			      std::vector<int> per_atoms_index = getBonded( molecules[i] , molecules[i].atom_index[j] );
			      int index_center = molecules[i].atom_index[j];
			      std::vector<double> position_center = getMinImage( atom_list , cell , index_center , index_center );
			      std::vector<double> position_atom1  = getMinImage( atom_list , cell , index_center , per_atoms_index[0] );
			      std::vector<double> position_atom2  = getMinImage( atom_list , cell , index_center , per_atoms_index[1] );
			      std::vector<double> position_atom3  = getMinImage( atom_list , cell , index_center , per_atoms_index[2] );
			      std::vector<double> vector_plan1 = difference( position_atom1 , position_atom2 ) ;
			      std::vector<double> vector_plan2 = difference( position_atom1 , position_atom3 ) ;
			      double dist = getDistanceFromPlan( vector_plan1 , vector_plan2 , position_center , position_atom1 );
			      if ( dist < 1.75 )
				{
				  distance_from_plan_in << step << " " << dist << std::endl;
				  distPlanC3_in.push_back( dist );
				}
			      co3_in++;
			    }
			}
		      else if ( angles.size() == 6 )
			{
			  if ( molecules[j].names.size() == 5 )
			    {
			      co4_alone++;
			    }
			  else co4_in++;
			}
		    }
		  else
		    {
		      if ( angles.size() == 1 )
			{
			  if( molecules[j].names.size() == 3 ) o2_in++;
			  else o2_alone++;
			}
		      else if ( angles.size() == 3 )
			{
			  if ( molecules[j].names.size() == 4 ) o3_in++;
			  else o3_alone++;
			}
		    }
		}
	    }
	  //-----------------------------------------------------------------------------
	  
	  //-----------------------------
	  // Stocking sizes of molecules
	  //----------------------------------------------------------
	  co2_stock_in.push_back( co2_in );
	  co3_stock_in.push_back( co3_in );
	  co4_stock_in.push_back( co4_in );
	  co2_stock_in_time << step << " " << co2_in << std::endl;
	  co3_stock_in_time << step << " " << co3_in << std::endl;
	  co4_stock_in_time << step << " " << co4_in << std::endl;
	  //------------------------------
	  co2_stock_alone.push_back( co2_alone );
	  co3_stock_alone.push_back( co3_alone );
	  co4_stock_alone.push_back( co4_alone );
	  co2_stock_alone_time << step << " " << co2_alone << std::endl;
	  co3_stock_alone_time << step << " " << co3_alone << std::endl;
	  co4_stock_alone_time << step << " " << co4_alone << std::endl;
	  //------------------------------
	  int c2_total = c2_alone + c2_in;
	  int c3_total = c3_alone + c3_in;
	  int c4_total = c4_alone + c4_in;
	  //------------------------------
	  std::cout << "c2_total: " << c2_total << std::endl;
	  std::cout << "c3_total: " << c3_total << std::endl;
	  std::cout << "c4_total: " << c4_total << std::endl;
	  double ratio=0;
	  if ( c2_total != 0 ) c2_ratio.push_back( (double)(c2_alone)/(double)(c2_total) );
	  else c2_ratio.push_back( ratio );
	  if ( c3_total != 0 ) c3_ratio.push_back( (double)(c3_alone)/(double)(c3_total) );
	  else c3_ratio.push_back( ratio );
	  if ( c4_total != 0 ) c4_ratio.push_back( (double)(c4_alone)/(double)(c4_total) );
	  else c4_ratio.push_back( ratio );
	  //-----------------------------------------------------------
	  o2_stock_in.push_back( o2_in );
	  o3_stock_in.push_back( o2_in );
	  o2_stock_alone.push_back( o2_alone );
	  o3_stock_alone.push_back( o3_alone );
	  //------------------------------
	  int o2_total = o2_alone + o2_in;
	  int o3_total = o3_alone + o3_in;
	  //------------------------------
	  if ( o2_total != 0 ) o2_ratio.push_back( (double)(o2_alone)/(double)(o2_total) );
	  else o2_ratio.push_back( ratio );
	  if ( o3_total != 0 ) o3_ratio.push_back( (double)(o3_alone)/(double)(o3_total) );
	  else o3_ratio.push_back( ratio );
	  //-----------------------------------------------------------

	  // Step making
	  //----------------------------------------------------
	  std::cout << step << std::endl;
	  //----------------------------------------------------
	}      
      step++;
     }
  //----------------------------------------------------

  //------------------
  // % co2 Histogram 
  //-------------------------------------------------------
  writeHistogram( co2_stock_in_out , normalizeHistogram( makeRegularHistogram( co2_stock_in , hist_co2_start, hist_co2_end, nb_co2_box ) ) );
  writeHistogram( co3_stock_in_out , normalizeHistogram( makeRegularHistogram( co3_stock_in , hist_co2_start, hist_co2_end, nb_co2_box ) ) );
  writeHistogram( co4_stock_in_out , normalizeHistogram( makeRegularHistogram( co4_stock_in , hist_co2_start, hist_co2_end, nb_co2_box ) ) );
  writeHistogram( co2_stock_alone_out , normalizeHistogram( makeRegularHistogram( co2_stock_alone , hist_co2_start, hist_co2_end, nb_co2_box ) ) );
  writeHistogram( co3_stock_alone_out , normalizeHistogram( makeRegularHistogram( co3_stock_alone , hist_co2_start, hist_co2_end, nb_co2_box ) ) );
  writeHistogram( co4_stock_alone_out , normalizeHistogram( makeRegularHistogram( co4_stock_alone , hist_co2_start, hist_co2_end, nb_co2_box ) ) );
  writeHistogram( c2_ratio_out , normalizeHistogram( makeRegularHistogram( c2_ratio , hist_start_ratio, hist_end_ratio, nb_box_ratio ) ) );
  writeHistogram( c3_ratio_out , normalizeHistogram( makeRegularHistogram( c3_ratio , hist_start_ratio , hist_end_ratio, nb_box_ratio ) ) );
  writeHistogram( c4_ratio_out , normalizeHistogram( makeRegularHistogram( c4_ratio , hist_start_ratio , hist_end_ratio, nb_co2_box ) ) );
  writeHistogram( distance_from_plan_hist_alone , normalizeHistogram( makeRegularHistogram( distPlanC3_alone , 0 , 1.75, 200 ) ) );
  writeHistogram( distance_from_plan_hist_in , normalizeHistogram( makeRegularHistogram( distPlanC3_in , -0.01 , 1.75, 200 ) ) );
  //-------------------------------------------------------
  
  //-----------------
  // Size histograms
  //---------------------------------------------------------
  for ( int i=0 ; i < hist.size() ; i++ )
    {
      hist[i].value = hist[i].value*center(hist[i]);
    }
  writeHistogram( hist_molecules , normalizeHistogram( makeRegularHistogram( sizes , hist_start , hist_end , nb_box ) ) );
  //----------------------------------------------------

  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  graph_molecules.close();
  hist_molecules.close();
  co2_stock_in_out.close();
  co3_stock_in_out.close();
  co4_stock_in_out.close();
  co2_stock_alone_out.close();
  co3_stock_alone_out.close();
  co4_stock_alone_out.close();
  c2_ratio_out.close();
  c3_ratio_out.close();
  c4_ratio_out.close();
  distance_from_plan_in.close();
  distance_from_plan_alone.close();
  distance_from_plan_hist_in.close();
  distance_from_plan_hist_alone.close();
  //----------------------
  
  return 0;
}
