**Boosting AIRSS**

Project led by Gabriele Mogni, supervised by Fabio Pietrucci

Aimed at studying the Ab Initio Random Search of Structure algorithm and to somehow create a faster method by making use of the PIV topological metric.

**Codes:**

*-castep.sh :* example of job to run castep on OCCIGEN (CINES)

*-airss_ex.sh :* example of job to run AIRSS on OCCIGEN (CINES)

*-treating_airss.sh :* short bash script that treats the result from an AIRSS launch and returns the energy and volume of the found structures

*-dpc_analysis :* script making use of the distances between selected structure to determine the convergence of the exploration of the configuration space ( written in Julia - requires some function in GPlib/Julia )