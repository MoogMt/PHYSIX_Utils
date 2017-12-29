bohr=0.5291772109;
hartree=4.35974417e-18;
ev=hartree/27.2;

function convertRy2eV{T <: Float64}(energy::T)
  return energy=energy*13.605693009
end
