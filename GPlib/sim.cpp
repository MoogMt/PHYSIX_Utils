#include "sim.h"

// Basic Sim
//================================================================================================
Sim emptySim() {
  std::vector<Atom> atoms;
 Cell cell = { 0, 0, 0 , 0, 0, 0 };
 Sim sim = { atoms , cell};
 return sim;
}
//================================================================================================

//===============
// Compress Box
//================================================================================================
Sim compressBox( Sim sim_set , double frac_a , double frac_b , double frac_c )
{
  sim_set.cell = compressBox( sim_set.cell , frac_a , frac_b , frac_c );
  sim_set.atoms = compressAtoms( sim_set.atoms, frac_a , frac_b , frac_c );
  return sim_set;
}
//================================================================================================


//================
// Diffusion Coef
//================================================================================================
void computeDiff( std::ofstream & output ,  std::ifstream & input , int comp_step , int start_step , int end_step , AtomList atom_list , AllTypeLUT lut_list , Cell cell )
{
  std::vector<double> x0,y0,z0;
  int step = 0;
  int count = 0;
  double d_coef = 0;
  AtomList atom_list0;
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step == start_step )
	{
	  x0 = atom_list.x;
	  y0 = atom_list.y;
	  z0 = atom_list.z;
	  atom_list0 = atom_list;
	}
      else if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  for ( int i=0 ; i < atom_list.x.size() ; i++ )
	    {
	      double dx = atom_list.x[i] - atom_list0.x[i];
	      double dy = atom_list.y[i] - atom_list0.y[i];
	      double dz = atom_list.z[i] - atom_list0.z[i];
	      dx = dx - (int)(dx);
	      dy = dy - (int)(dy);
	      dz = dz - (int)(dz);
	      atom_list.x[i] = atom_list0.x[i] + dx;
	      atom_list.y[i] = atom_list0.y[i] + dy;
	      atom_list.z[i] = atom_list0.z[i] + dz;
	      atom_list0.x[i] = atom_list.x[i];
	      atom_list0.y[i] = atom_list.y[i];
	      atom_list0.z[i] = atom_list.z[i];
	    }
	  std::vector<double> x = difference( atom_list.x , x0 ) ;
	  std::vector<double> y = difference( atom_list.y , y0 );
	  std::vector<double> z = difference( atom_list.z , z0 );
	  std::vector<double> r; r.assign( x.size() , 0 );
	  for ( int i=0; i < x.size() ; i++ )
	    {
	      r[i] = abs(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
	    }
	  output << (step-start_step) << " " << average(r) << std::endl;
	}      
      std::cout << step << std::endl;
      step++;
    }
  return;
}
//================================================================================================
