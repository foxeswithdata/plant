// -*-c++-*-
#ifndef PLANT_PLANT_NEW_STRATEGY_H_
#define PLANT_PLANT_NEW_STRATEGY_H_

#include <memory>
#include <plant/control.h>
#include <plant/qag_internals.h> // quadrature::intervals_type
#include <plant/internals.h> // quadrature::intervals_type
#include <plant/strategy.h>

#include <plant/models/NEW_environment.h>1

namespace plant {

class NEW_Strategy: public Strategy<NEW_Environment> {
public:
  typedef std::shared_ptr<NEW_Strategy> ptr;
  NEW_Strategy(); // Constructor definition

  /***
   * update this when the length of state_names changes
   * replace with the number of states that the strategy tracks
   */
  static size_t state_size () { return 3; }

  /***
   * update this when the length of aux_names changes
   */
  size_t aux_size () { return aux_names().size(); }

  /***
   * Replace with relevant states in the plant strategy
   */
  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "state_1",
      "state_2",
      "state_3"
    });
  }

  /***
   * don't understand this yet
   */
  std::vector<std::string> aux_names() {
    std::vector<std::string> ret({
      "competition_effect",
      "net_mass_production_dt"
    });
    // add the associated computation to compute_rates and compute there
    if (collect_all_auxillary) {
      ret.push_back("area_sapwood");
    }
    return ret;
  }


  /***
   * Create mapping between state and aux names and indices
   * in state and aux files
   */
  void refresh_indices();

  /***
   * This method is called to compute the state ODEs
   * This is where all the behaviour method will be called
   * Can the method arguments be changed?
   */
  void compute_rates(const NEW_Environment& environment, bool reuse_intervals,
                     Internals& vars);

  /***
   * still don't understand this
   */
  void update_dependent_aux(const int index, Internals& vars);


  /***
   * Define all the behaviours that will be met by the stratgegy
   * Why const?
   */
  double equation_1(double arg_1, arg_2);

  /***
   * Set constants within the strategy
   */
  void prepare_strategy();


  /***
   * Define parameters etc.
   */
  double par1, par2, par3;

  /***
   * Name of the strategy
   */
  std::string name;

  /***
   * Still need to figure this out
   */
  Assimilation assimilator;

};


/***
 * Make a strategy pointer
 */
NEW_Strategy::ptr make_strategy_ptr(NEW_Strategy s);

}

#endif
