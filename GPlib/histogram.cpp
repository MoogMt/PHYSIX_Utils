#include "histogram.h"

// BINS
//-------------------------------------------------
Bin emptyBin()
{
  Bin bin = { 0, 0, 0};

  return bin;
}
//
Bin makeBinReal( double begin , double end , int value )
{
  Bin bin = { begin , end , value };
  return bin;
}
//
double center( Bin bin )
{
  return (bin.end+bin.begin)/2.0;
}
//
bool overlap( Bin bin1, Bin bin2 )
{
  if ( max(bin1.begin,bin2.begin) <  min(bin1.end,bin2.end) )
    {
      return true;
    }
  else
    {
      return false;
    }
}
//
void fillBin( Bin &bin , std::vector<double> data )
{
  for ( int i=0 ; i < data.size() ; i++ )
    {
      if ( data[i] > bin.begin && data[i] < bin.end )
	{
	  bin.value++;
	}
    }
  return ;
}
//
Bin makeBin( double bin_min, double bin_max, std::vector<double> data)
{
  Bin bin = { bin_min , bin_max , 0 };
  fillBin( bin , data );
  return bin;
}
//
Bin addBinsMin( Bin bin1, Bin bin2)
{
  Bin bin_final;
  if ( overlap(bin1,bin2) )
    {
      bin_final.begin = max(bin1.begin,bin2.begin);
      bin_final.end   = min(bin1.end,bin2.end);
      bin_final.value = bin1.value + bin2.value;
      return bin_final;
    }
  else
    {
      return emptyBin();
    }
}
//
Bin addBinsMax( Bin bin1, Bin bin2 )
{
  Bin bin_final = { min(bin1.begin,bin2.begin) , max(bin1.end,bin2.end) , bin1.value + bin2.value };
  return bin_final;
}
//-------------------------------------------------

// MAKE HISTOGRAM
//--------------------------------------------------------------------------------------------------
std::vector<Bin> makeRegularHistogram( std::vector<double> data , double x_min , double x_max , int number_bins )
{
  std::vector<Bin> bins_hist;
  double delta = (x_max-x_min)/(double)(number_bins);
  if ( delta < 0 )
    {
      std::cout << "Error: Min of the interval to max!" << std::endl;
      return bins_hist;
    }
  for ( int i=0 ; i < number_bins ; i++ )
    {
      bins_hist.push_back( makeBin( x_min + i*delta , x_min + (i+1)*delta , data ) );
    }
  return bins_hist;
}
//
std::vector<Bin> makeRegularHistogram( std::vector<double> data_x , std::vector<double> data_y , int number_bins )
{
  std::vector<Bin> bins_hist;
  double x_min = min(data_x);
  double x_max = max(data_x);
  double delta = (x_max-x_min)/(double)(number_bins);
  for ( int i=0 ; i < number_bins ; i++ )
    {
      bins_hist.push_back( makeBin( x_min + i*delta , x_min + (i+1)*delta , data_y ) );
    }
  return bins_hist;
}
//
std::vector<Bin> makeHistograms( std::vector<double> data, std::vector<double> bins_limits)
{
  std::vector<Bin> bins_hist;
  for( int i=0 ; i < bins_limits.size()-1 ; i++ )
    {
      bins_hist.push_back( makeBin( bins_limits[i] , bins_limits[i+1] , data ) );
    }
  return bins_hist;
}
std::vector<Bin> makeHistograms( std::vector<double> data, std::vector<Bin> bins)
{
  for( int i=0 ; i < bins.size() ; i++ )
    {
      fillBin( bins[i] , data ) ;
    }
  return bins;
}
//---------------------------------------------------------------------------------------
// Modification
//---------------------------------------------------------------------------------------
std::vector<Bin> addHistograms( std::vector<Bin> hist1 , std::vector<Bin> hist2 )
{
  std::vector<Bin> sum_hist;
  for ( int i=0 ; i < hist1.size() ; i++ )
    {
      for ( int j=0 ; j  < hist2.size() ; j++ )
	{
	  if ( overlap( hist1[i] , hist2[j] ) )
	    {
	      sum_hist.push_back( addBinsMin( hist1[i] , hist2[j] ) );
	    }
	}
    }
  return sum_hist;
}
//---------------------------------------------------------------------------------------
// Total Value
//---------------------------------------------------------------------------------------
double getTotalValue( std::vector<Bin> hist )
{
  double total_value = 0;
  for ( int i=0 ; i < hist.size() ; i++ )
    {
      total_value += hist[i].value;
    }
  return total_value;
}
//-----
// IO
//---------------------------------------------------------------------------------------
void writeHistogram( std::ofstream & file , std::vector<Bin> hist )
{
  for ( int i = 0 ; i < hist.size() ; i++ )
    {
      file  << center(hist[i]) << " " << hist[i].value << std::endl;
    }
  return;
}
void writeHistogram( std::string file_name , std::vector<Bin> hist )
{
  std::ofstream file ( file_name.c_str() ,  std::ios::out | std::ios::app );
  for ( int i=0 ; i < hist.size() ; i++ )
    {
      file << center(hist[i]) << " " << hist[i].value << std::endl;
    }
  file.close();
  return;
}
//---------------------------------------------------------------------------------------
//==========================================================================================


// Bin Real
//-------------------------------------------------------------------------------------------------
BinReal emptyBinReal()
{
  BinReal bin = { 0, 0, 0};
  return bin;
}
//
BinReal makeBinReal( double begin , double end , double value )
{
  BinReal bin = { begin , end , value };
  return bin;
}
//
double center( BinReal bin )
{
  return ( bin.begin + bin.end )/2;
}
//
bool overlap( BinReal bin1 , BinReal bin2 )
{ 
  if ( max(bin1.begin,bin2.begin) <  min(bin1.end,bin2.end) )
    {
      return true;
    }
  else
    {
      return false;
    }
}
//
void fillBin( BinReal &bin , std::vector<double> data )
{
  for ( int i=0 ; i < data.size() ; i++ )
    {
      if ( data[i] > bin.begin && data[i] < bin.end )
	{
	  bin.value++;
	}
    }
  return ;
}
//
BinReal makeBinReal( double bin_min, double bin_max, std::vector<double> data)
{
  BinReal bin = { bin_min , bin_max , 0 };
  fillBin( bin , data );
  return bin;
}
//
BinReal addBinsMin( BinReal bin1, BinReal bin2)
{
  BinReal bin_final;
  if ( overlap(bin1,bin2) )
    {
      bin_final.begin = max(bin1.begin,bin2.begin);
      bin_final.end   = min(bin1.end,bin2.end);
      bin_final.value = bin1.value + bin2.value;
      return bin_final;
    }
  else
    {
      return emptyBinReal();
    }
}
//
BinReal addBinsMax( BinReal bin1, BinReal bin2 )
{
  BinReal bin_final = { min(bin1.begin,bin2.begin) , max(bin1.end,bin2.end) , bin1.value + bin2.value };
  return bin_final;
}
//-------------------------------------------------------------------------------------------------

// Bin Real Histogram
//-------------------------------------------------------------------------------------------------
std::vector<BinReal> normalizeHistogram( std::vector<Bin> hist )
{
  std::vector<BinReal> hist_real;
  double norm = getTotalValue( hist );
  for ( int i=0 ; i < hist.size() ; i++ )
    {
      hist_real.push_back( makeBinReal( hist[i].begin , hist[i].end , hist[i].value/norm ) );
    }
  return hist_real;
}
//
void writeHistogram( std::ofstream & file , std::vector<BinReal> hist )
{
  for ( int i=0 ; i < hist.size() ; i++ )
    {
      file << center(hist[i]) << " " << hist[i].value << std::endl;
    }
  return;
}
//
void writeHistogram( std::string file_name , std::vector<BinReal> hist )
{
  std::ofstream file ( file_name.c_str() ,  std::ios::out | std::ios::app );
  for ( int i=0 ; i < hist.size() ; i++ )
    {
      file << center(hist[i]) << " " << hist[i].value << std::endl;
    }
  file.close();
  return;
}
//
bool checkSizeHists( std::vector<BinReal> hist , std::vector<BinReal> hist2 )
{
  if ( hist.size() == hist2.size() )
    {
      return true;
    }
  else
    {
      return false;
    }
}
//
bool checkSizeHists( std::vector< std::vector<BinReal> > hist_list )
{
  for ( int i=0 ; i < hist_list.size()-1 ; i++ )
    {
      for ( int j=0 ; j < hist_list.size() ; j++ )
	{
	  if ( ! checkSizeHists( hist_list[i] , hist_list[j] ) )
	    {
	      return false;
	    }
	}
    }
  return true;
}
//
void writeBinReal( std::ofstream & file , BinReal bin, bool wcenter)
{
  if ( wcenter )
    {
      file << center(bin) <<  " " << bin.value << " ";
    }
  else
    {
      file << bin.value << " ";
    }
  return ;
}
//
void writeBinRealCenter( std::ofstream & file, BinReal bin )
{
  file << center(bin) << " ";
  return;
}
//
void writeHistBinCenter( std::ofstream & file , std::vector<BinReal> hist , int index )
{
  writeBinRealCenter( file, hist[index] );
}  
//
void writeHistBin( std::ofstream & file , std::vector<BinReal> hist, int index , bool wcenter )
{
  writeBinReal( file , hist[index] , wcenter );
}
//
void writeHistograms( std::ofstream & file , std::vector< std::vector<BinReal> > hist_list )
{
  if ( checkSizeHists( hist_list ) )
    {
      std::cout << "Error: Some histograms do not have the same size." << std::endl;
      return;
    }
  for ( int i=0 ; i < hist_list[0].size() ; i++ )
    {
      writeHistBinCenter( file, hist_list[0] , i );
      for ( int j=0 ; j < hist_list.size() ; j++ )
	{
	  writeHistBin( file , hist_list[j] , i , false );
	}
      file << std::endl;
    }  
  return ;
}
//-------------------------------------------------------------------------------------------------



