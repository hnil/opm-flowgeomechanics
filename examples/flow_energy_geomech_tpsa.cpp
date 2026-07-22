/*
  Fracture-enabled thermal water injection simulator with the TPSA
  geomechanics backend.  Companion to flow_energy_geomech (VEM backend);
  requires MECHSOLV in the deck to select the TPSA solver and coupling
  scheme.
*/
#include "config.h"

#include <opm/simulators/flow/Main.hpp>

#include "MechTpsaTypeTag.hpp"

int
main(int argc, char** argv)
{
    OPM_TIMEBLOCK(fullSimulation);
    using TypeTag = Opm::Properties::TTag::FlowProblemMechTpsa;
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
}
