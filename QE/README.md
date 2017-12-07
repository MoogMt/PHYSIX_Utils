STUFF ABOUT INSTALLING QE:

==========================
-> Weird Compilation Stuff
===========================

- ALWAYS use the method "more is better" when loading modules. The more modules you load, the best chance you have to
see the compilation suceed.

- If trouble when installing simple pw:
make veryclean 
./configure --with-scalpack="intel"

Sometimes helps, as it uses scalapack from intel not mkl which tends to makes the compilation go banana

- NEVER use the -j N option when doing make plumed, this will make you think it does patch, when actually it does not...
(lost lots of time due to this qwerk)
