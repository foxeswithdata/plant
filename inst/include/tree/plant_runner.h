// -*-c++-*-
#ifndef TREE_TREE_PLANT_RUNNER_H_
#define TREE_TREE_PLANT_RUNNER_H_

#include <tree/ffw16_plant_plus.h>
#include <tree/environment.h>

namespace tree {
namespace tools {

// TODO: This should be templated, I think, but that plays badly with
// RcppR6's requirements for "concrete" types.
struct PlantRunner {
  PlantRunner(FFW16_PlantPlus plant_, Environment environment_)
    : plant(plant_), environment(environment_) {
    plant.compute_vars_phys(environment);
  }
  static size_t ode_size() {return FFW16_PlantPlus::ode_size();}
  double ode_time() const {return environment.time;}
  ode::const_iterator set_ode_state(ode::const_iterator it, double time) {
    it = plant.set_ode_state(it);
    environment.time = time;
    plant.compute_vars_phys(environment);
    return it;
  }
  ode::iterator ode_state(ode::iterator it) const {
    return plant.ode_state(it);
  }
  ode::iterator ode_rates(ode::iterator it) const {
    return plant.ode_rates(it);
  }
  FFW16_PlantPlus plant;
  Environment environment;
};

}
}

#endif
