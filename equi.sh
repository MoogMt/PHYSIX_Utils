#!/bin/csh

# E=f(P) batcher


# Parameters
ecut = "80"     # Cut-Off energy for wavefunction
ecutrho = "500" # Cut-Off energy for density

while read line
    press = $(echo $line | cut -d " " -f1 
    set coorda = 5.614332 # a en u.a.
    set coordb = 1.       # b/a 
    set coordc = 1.042    # c/a
    
    set -x
    
    pwd=$(pwd)
	
    cp ./* $TMPDIR
    
	######################### VC-RELAX ##############################
    
	cat > pw.d << ***
	
    &control
    calculation = 'vc-relax',
    restart_mode='from_scratch',
    prefix='amf.$run',
	pseudo_dir='./',
  tstress = .true.
	tprnfor = .true.
    outdir='./'
    verbosity = 'high'
    /
    &system
    ibrav=0, 
	celldm(1)=$coorda,
    celldm(2)=$coordb,
    celldm(3)=$coordc,
    nat=6, ntyp=2,
    ecutwfc=$ecut,
    ecutrho=$ecutrho,
	/
    &electrons
    mixing_mode = 'plain'
    mixing_beta = 0.3
    conv_thr =  1.0d-12
    /
    &ions
    /
    &cell                
    press=$press
    /   
	ATOMIC_SPECIES
	C 1.0079  C.pbe-van_ak.UPF
	O 14.0067 O.pbe-van_ak.UPF
    ATOMIC_POSITIONS crystal
	C        0.000000000   0.000000000   0.000000000
    O        0.500000000   0.500000000   0.500000000
    O        0.316070721   0.316070721   0.316070721
    C        1.316070721   1.683929279   0.68392927
    O        0.683929279   1.316070721   1.683929279
    O        1.683929279   0.683929279   1.316070721
	K_POINTS automatic
    6 6 4 1 1 1
    ***
    
    poe pw.x < pw.d > NH4F_P-42m.relax.$ecut.$ecutrho.$press.out
	    
done < PresList

