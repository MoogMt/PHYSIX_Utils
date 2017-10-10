#include "sim.h"

Sim emptySim() {
  std::vector<Atom> atoms;
 Cell cell = { 0, 0, 0 , 0, 0, 0 };
 Sim sim = { atoms , cell};
 return sim;
}
