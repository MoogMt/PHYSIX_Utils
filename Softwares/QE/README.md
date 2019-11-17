STUFF ABOUT INSTALLING QE:

===Weird Compilation Stuff===

- NEVER use the -j N option when doing make plumed, this will make you think it does patch, when actually it does not...
(lost lots of time due to this qwerk)

- PLUMED does not work with QE for version 6.x onwards , revert to 5.x 

- ALWAYS use the method "more is better" when loading modules. The more modules you load, the best chance you have to
see the compilation suceed.

- If trouble when installing simple pw:
make veryclean 
./configure --with-scalpack="intel"

- On MeSU, it is necessary to add an intel lib for PLUMED or compilation will crash

- On Ada:

  -> Always use the --with-scalapack="intel" option for ./configure or compilation will crash
  
  -> Add latest gcc compiler module for PLUMED or compilation will crash


