#include "histogram.h"

// BINS
//-------------------------------------------------
void fillBin( Bin & bin , std::vector<double> data )
{
  for ( int i=0 ; i < data.size() ; i++ )
    {
      if ( data[i] > min && data[i] < max )
	{
	  bin.value++;
	  i--;
	  data.erase(erase.begin() + i );
	}
    }
  return ;
}
Bin makeBin( int bin_min, int bin_max, std::vector<double> & data)
{
  Bin bin = { bin_min , bin_max , 0 };
  return fillBin( bin& , data& );
}
//-------------------------------------------------

// MAKE HISTOGRAM
//--------------------------------------------------------------------------------------------------
std::vector<bins> makeRegularHistogram( std::vector<double> data_x , std::vector<double> data_y , double x_min , double x_max , int number_bins )
{
  std::vector<Bin> bins_hist;
  for ( int i=0 ; i < number_bins ; i++ )
    {
      bins_hist.push_back( makeBin( ) );
    }
  returns bins_hist;
}
//--------------------------------------------------------------------------------------------------
