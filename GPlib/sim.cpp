#include "sim.h"

Sim emptySim() {
  std::vector<Atom> atoms;
 Cell cell = { 0, 0, 0 , 0, 0, 0 };
 Sim sim = { atoms , cell};
 return sim;
}

Sim compressBox( Sim sim_set , double frac_a , double frac_b , double frac_c )
{
  sim_set.cell = compressBox( sim_set.cell , frac_a , frac_b , frac_c );
  sim_set.atoms = compressAtoms( sim_set.atoms, frac_a , frac_b , frac_c );
  return sim_set;
}

void computeDiffCoef( std::ofstream & output ,  std::ifstream & input , int comp_step , int start_step , int end_step , AtomList atom_list , AllTypeLUT lut_list , Cell cell )
{
  std::vector<double> x0,y0,z0;
  int step = 0;
  int count=0;
  double d_coef = 0;
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step == start_step )
	{
	  x0 = backIn( atom_list.x , cell.a );
	  y0 = backIn( atom_list.y , cell.b );
	  z0 = backIn( atom_list.z , cell.c );
	}
      else if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  std::vector<double> x = difference( backIn( atom_list.x, cell.a ) , x0 );
	  std::vector<double> y = difference( backIn( atom_list.y, cell.b ) , y0 );
	  std::vector<double> z = difference( backIn( atom_list.z, cell.c ) , z0 );
	  std::vector<double> r = square( squaroot( addVector( addVector( square( x ), square( y ) ), square( z ) ) ) );
	  output << step-start_step << " " << average( r ) << std::endl;
	  count++;
	  if ( step > start_step + 2000 )
	    {
	      d_coef += average( r );
	    }
	}      
      std::cout << step << std::endl;
      step++;
    }
  return;
}
