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
  std::ifstream input("histogram.dat");
  //--------------------------------

  //--------
  // Output
  //-------------------------------------
  std::ofstream output2("angles2.dat");
  std::ofstream output3("angles3.dat");
  std::ofstream output4("angles4.dat");
  //-------------------------------------
  
  //-----------------
  // Read histograms
  //---------------------------------------------------------------------------------------------
  // CO2
  //---------------------------------------------------------------------------------------------
  std::vector<BinReal> angles2_882_2000_hist = readRegularHistogramReal("8.82/2000K/pressure.dat");
  std::vector<BinReal> angles2_900_2000_hist = readRegularHistogramReal("9.0/2000K/pressure.dat");
  std::vector<BinReal> angles2_980_2000_hist = readRegularHistogramReal("9.8/2000K/pressure.dat");
  std::vector<BinReal> angles2_882_2500_hist = readRegularHistogramReal("8.82/2500K/pressure.dat");
  std::vector<BinReal> angles2_900_2500_hist = readRegularHistogramReal("9.0/2500K/pressure.dat");
  std::vector<BinReal> angles2_980_2500_hist = readRegularHistogramReal("9.8/2500K/pressure.dat");
  std::vector<BinReal> angles2_882_3000_hist = readRegularHistogramReal("8.82/3000K/pressure.dat");
  std::vector<BinReal> angles2_900_3000_hist = readRegularHistogramReal("9.0/3000K/pressure.dat");
  std::vector<BinReal> angles2_980_3000_hist = readRegularHistogramReal("9.8/3000K/pressure.dat");
  //---------------------------------------------------------------------------------------------
  // CO3
  //---------------------------------------------------------------------------------------------
  std::vector<BinReal> angles3_882_2000_hist = readRegularHistogramReal("8.82/2000K/pressure.dat");
  std::vector<BinReal> angles3_900_2000_hist = readRegularHistogramReal("9.0/2000K/pressure.dat");
  std::vector<BinReal> angles3_980_2000_hist = readRegularHistogramReal("9.8/2000K/pressure.dat");
  std::vector<BinReal> angles3_882_2500_hist = readRegularHistogramReal("8.82/2500K/pressure.dat");
  std::vector<BinReal> angles3_900_2500_hist = readRegularHistogramReal("9.0/2500K/pressure.dat");
  std::vector<BinReal> angles3_980_2500_hist = readRegularHistogramReal("9.8/2500K/pressure.dat");
  std::vector<BinReal> angles3_882_3000_hist = readRegularHistogramReal("8.82/3000K/pressure.dat");
  std::vector<BinReal> angles3_900_3000_hist = readRegularHistogramReal("9.0/3000K/pressure.dat");
  std::vector<BinReal> angles3_980_3000_hist = readRegularHistogramReal("9.8/3000K/pressure.dat");
  //---------------------------------------------------------------------------------------------
  // CO4
  //---------------------------------------------------------------------------------------------
  std::vector<BinReal> angles4_882_2000_hist = readRegularHistogramReal("8.82/2000K/pressure.dat");
  std::vector<BinReal> angles4_900_2000_hist = readRegularHistogramReal("9.0/2000K/pressure.dat");
  std::vector<BinReal> angles4_980_2000_hist = readRegularHistogramReal("9.8/2000K/pressure.dat");
  std::vector<BinReal> angles4_882_2500_hist = readRegularHistogramReal("8.82/2500K/pressure.dat");
  std::vector<BinReal> angles4_900_2500_hist = readRegularHistogramReal("9.0/2500K/pressure.dat");
  std::vector<BinReal> angles4_980_2500_hist = readRegularHistogramReal("9.8/2500K/pressure.dat");
  std::vector<BinReal> angles4_882_3000_hist = readRegularHistogramReal("8.82/3000K/pressure.dat");
  std::vector<BinReal> angles4_900_3000_hist = readRegularHistogramReal("9.0/3000K/pressure.dat");
  std::vector<BinReal> angles4_980_3000_hist = readRegularHistogramReal("9.8/3000K/pressure.dat");
  //---------------------------------------------------------------------------------------------

  //------------------
  // Compute averages
  //---------------------
  std::cout << "8.82 " << average(press_882_2000_hist) << " " << average(press_882_2500_hist) << " " << average(press_882_3000_hist) << " " << std::endl;
  std::cout << "9.00 " << average(press_900_2000_hist) << " " << average(press_900_2500_hist) << " " << average(press_900_3000_hist) << " " << std::endl;
  std::cout << "9.80 " << average(press_980_2000_hist) << " " << average(press_980_2500_hist) << " " << average(press_980_3000_hist) << " " << std::endl;
  //---------------------
    
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  output.close();
  //----------------------
  
  return 0;
}
