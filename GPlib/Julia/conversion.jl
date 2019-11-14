module conversion

AU2Pascal=2.9421912e13
Pascal2AU=1/AU2Pascal

THz2cm=0.0299792458
cm2THz=1/0.0299792458

Hz2cm=33.3565
cm2Hz=1/33.3565

hartree=4.35974417e-18;

function ry2eV( energy::T ) where { T <: Real }
  return energy*13.605693009
end

function ev2Ry( energy::T ) where { T <: Real }
  return energy/13.605693009
end

function ha2eV( energy::T ) where { T <: Real }
  return energy/27.2
end

function bohr2Ang( distance::T ) where { T <: Real }
  return distance*0.5291772109
end

function Ang2Bohr( distance::T ) where { T <: Real }
  return distance/0.5291772109
end

end
