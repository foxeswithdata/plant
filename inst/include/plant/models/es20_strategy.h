// -*-c++-*-
#ifndef PLANT_PLANT_ES20_STRATEGY_H_
#define PLANT_PLANT_ES20_STRATEGY_H_

#include <memory>
#include <plant/control.h>
#include <plant/qag_internals.h> // quadrature::intervals_type
#include <plant/internals.h> // quadrature::intervals_type
#include <plant/strategy.h>
#include <plant/models/ff16_strategy.h> // extending this strategy
#include <plant/models/ff16_environment.h>
#include <plant/models/assimilation.h>

namespace plant {

class ES20_Strategy: public FF16_Strategy, public Strategy<ES20_Environment> {
public:
  typedef std::shared_ptr<ES20_Strategy> ptr;
  ES20_Strategy();

  // update this when the length of state_names changes
//  static size_t state_size () { return 5; }
  // update this when the length of aux_names changes
//  size_t aux_size () { return aux_names().size(); }

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
//      "mortality",
//      "fecundity",
      "area_heartwood",
      "mass_heartwood",

      "mass_storage",

      "area_leaf",
      "mass_leaf",

      "area_sapwood",
      "mass_sapwood",

      "area_bark",
      "mass_bark",

      "mass_root",

     "area_stem",

      "diameter"
      });
  }


ES20_Strategy::ptr make_strategy_ptr(ES20_Strategy s);

}

#endif

