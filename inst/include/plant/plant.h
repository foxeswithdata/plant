// -*-c++-*-
#ifndef PLANT_PLANT_PLANT_MINIMAL_H_
#define PLANT_PLANT_PLANT_MINIMAL_H_

#include <memory> // std::shared_ptr
#include <plant/ode_interface.h>
#include <vector>
#include <plant/internals.h>
#include <plant/util.h>
#include <plant/uniroot.h>

#include <plant/models/ff16_strategy.h>
#include <plant/models/ff16r_strategy.h>
#include <plant/models/es20_strategy.h>
#include <plant/models/es20r_strategy.h>

#include <iostream>


namespace plant {

template <typename T, typename E> class Plant {
public:
  typedef T strategy_type;
  typedef E environment_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  // for the time being...
  Plant(strategy_type_ptr s) : strategy(s) {
    if (strategy->aux_index.size() != s->aux_size()) {
      strategy->refresh_indices();
    }
    vars.resize(strategy_type::state_size(), s->aux_size()); // = Internals(strategy_type::state_size());
    set_state("height", strategy->height_0);
  }
  
  // useage: state(HEIGHT_INDEX)
  double state(std::string name) const {
    return vars.state(strategy->state_index.at(name));
  }
  double state(int i) const { return vars.state(i); }
  
  // useage:_rate("area_heartwood")
  double rate(std::string name) const {
    return vars.rate(strategy->state_index.at(name));
  }
  double rate(int i) const { return vars.rate(i); }

  // useage: set_state("height", 2.0)
  void set_state(std::string name, double v) {
    int i = strategy->state_index.at(name);
    vars.set_state(i, v);
    strategy->update_dependent_aux(i, vars);
  }
  void set_state(int i, double v) {
    vars.set_state(i, v);
    strategy->update_dependent_aux(i, vars);
  }

  // aux vars by name and index
  double aux(std::string name) const {
    return vars.aux(strategy->aux_index.at(name));
  }
  double aux(int i) const { return vars.aux(i); } 

  double compute_competition(double z) const {
    return strategy->compute_competition(z, state(HEIGHT_INDEX)); // aux("competition_effect"));
  }

  void compute_rates(const environment_type &environment,
                         bool reuse_intervals = false) {
    strategy->compute_rates(environment, reuse_intervals, vars);
  }
  
  double establishment_probability(const environment_type &environment) {
    return strategy->establishment_probability(environment);
  }

  double net_mass_production_dt(const environment_type &environment) {
    // TODO:  maybe reuse intervals? default false 
    // should fix with a switch 

      if(strategy_name() == "ES20"){
        return strategy->net_mass_production_dt(environment, state(HEIGHT_INDEX), 
                                                     state(AREA_LEAF_INDEX), state(MASS_LEAF_INDEX),
                                                     state(MASS_SAPWOOD_INDEX), state(MASS_BARK_INDEX),
                                                     state(MASS_ROOT_INDEX), aux("competition_effect"));
      }
      else if(strategy_name() == "ESr20"){
        return strategy->net_mass_production_dt(environment, state(HEIGHT_INDEX), 
                                                state(AREA_LEAF_INDEX), state(MASS_LEAF_INDEX),
                                                state(MASS_SAPWOOD_INDEX), state(MASS_BARK_INDEX),
                                                state(MASS_ROOT_INDEX), aux("competition_effect"));
      }
      else if(strategy_name() == "FF16"){
        return strategy->net_mass_production_dt(environment, state(HEIGHT_INDEX), 
                                                state(AREA_LEAF_INDEX), strategy->mass_leaf(state(AREA_LEAF_INDEX)),
                                                strategy->mass_sapwood(strategy->area_sapwood(state(AREA_LEAF_INDEX)), state(HEIGHT_INDEX)), 
                                                strategy->mass_bark(strategy->area_bark(state(AREA_LEAF_INDEX)), state(HEIGHT_INDEX)),
                                                strategy->mass_root(state(AREA_LEAF_INDEX)),aux("competition_effect"));
      }
      else if(strategy_name() == "FF16r"){
        return strategy->net_mass_production_dt(environment, state(HEIGHT_INDEX), 
                                                state(AREA_LEAF_INDEX), strategy->mass_leaf(state(AREA_LEAF_INDEX)),
                                                strategy->mass_sapwood(strategy->area_sapwood(state(AREA_LEAF_INDEX)), state(HEIGHT_INDEX)), 
                                                strategy->mass_bark(strategy->area_bark(state(AREA_LEAF_INDEX)), state(HEIGHT_INDEX)),
                                                strategy->mass_root(state(AREA_LEAF_INDEX)),aux("competition_effect"));
      }
      else{
        return 0;
      }
   
  }

  // * ODE interface
  static size_t ode_size() { return strategy_type::state_size(); }
  static std::vector<std::string> ode_names() { return strategy_type::state_names(); }

  size_t aux_size() { return strategy->aux_size(); }
  std::vector<std::string> aux_names() { return strategy->aux_names(); }

  ode::const_iterator set_ode_state(ode::const_iterator it) {
    for (int i = 0; i < vars.state_size; i++) {
      vars.states[i] = *it++;
      strategy->update_dependent_aux(i, vars);
    }
    return it;
  }
  ode::iterator ode_state(ode::iterator it) const {
    for (int i = 0; i < vars.state_size; i++) {
      *it++ = vars.states[i];
    }
    return it;
  }
  ode::iterator ode_rates(ode::iterator it) const {
    for (int i = 0; i < vars.state_size; i++) {
      *it++ = vars.rates[i];
    }
    return it;
  }

  // Single individual methods

  // Used in the stochastic model:
  double mortality_probability() const { return 1 - exp(-state(MORTALITY_INDEX)); }
  
  void reset_mortality() { set_state("mortality", 0.0); }

  // TODO: Eventually change to growth rate given size
  double growth_rate_given_height(double height, const environment_type& environment) {
    set_state("height", height);
    compute_rates(environment, true);
    return rate("height");
  }

  // TODO recame lcp_whole_plant to competition_compensation_point
  double lcp_whole_plant() {
    environment_type env = environment_type();

    auto target = [&] (double x) mutable -> double {
      env.set_fixed_environment(x);
      compute_rates(env);
      return net_mass_production_dt(env);
    };

    const double f1 = target(1.0);
    if (f1 < 0.0) {
      return NA_REAL;
    } else {
      const double tol = control().plant_seed_tol;
      const size_t max_iterations = control().plant_seed_iterations;
      return util::uniroot(target, 0.0, 1.0, tol, max_iterations);
    }
  }



  std::string strategy_name() const { return strategy->name; }

  // * R interface
  strategy_type r_get_strategy() const { return *strategy.get(); }

  // ! External R code depends on knowing r internals for like growing plant to
  // ! height or something
  Internals r_internals() const { return vars; }
  const Control &control() const { return strategy->control; }
  
  // convert string to integer for use in switch statements (source: https://stackoverflow.com/questions/16388510/evaluate-a-string-with-a-switch-in-c)
  // move this to UTIL to make generally available
  constexpr unsigned int str2int(const char* str, int h = 0) const {
    return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
  }

private:
  strategy_type_ptr strategy;
  Internals vars;
};

template <typename T, typename E> Plant<T,E> make_plant(T s) {
  return Plant<T,E>(make_strategy_ptr(s));
}



} // namespace plant

#endif
