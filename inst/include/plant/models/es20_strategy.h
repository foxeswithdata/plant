// Built from  inst/include/plant/models/ff16_strategy.h on Wed Sep 16 07:11:11 2020 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_ES20_STRATEGY_H_
#define PLANT_PLANT_ES20_STRATEGY_H_

#include <memory>
#include <plant/strategy.h>
#include <plant/models/es20_environment.h>
#include <plant/models/es20_assimilation.h>

namespace plant {

class ES20_Strategy: public Strategy<ES20_Environment> {
public:
  typedef std::shared_ptr<ES20_Strategy> ptr;
  ES20_Strategy();

  // Overrides ----------------------------------------------

  // update this when the length of state_names changes
  static size_t state_size () { return 15; }
  // update this when the length of aux_names changes
  size_t aux_size () { return aux_names().size(); }

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
      "mortality",
      "fecundity",
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
      "diameter_stem"
      });
  }

  std::vector<std::string> aux_names() {
    std::vector<std::string> ret({
      "competition_effect",
      "net_mass_production_dt",
      "dbiomass_dt",
      "respiration_dt",
      "respiration_leaf_dt",
      "respiration_sapwood_dt",
      "respiration_root_dt",
      "respiration_bark_dt",
      "darea_leaf_dmass_live",
      "darea_leaf_dt",
      "dmass_leaf_darea_leaf",
      "dmass_sapwood_darea_leaf",
      "dmass_bark_darea_leaf",
      "dmass_root_darea_leaf",
      "dmass_leaf_dt",
      "dmass_sapwood_dt",
      "dmass_bark_dt",
      "dmass_root_dt",
      "dheight_darea_leaf",
      "height_difference",
      "height_adjustment", 
      "mass_sapwood_difference",
      "mass_sapwood_adjustment",
      "fecundity_dt_abs",
      "storage_portion",
      "dead_flag",
      "mortality_storage_dependent_dt",
      "mortality_growth_dependent_dt",
      "productivity_area",
      "stress"
    });
    return ret;
  }

  // Translate generic methods to ES20 strategy leaf area methods

  double competition_effect(double height) const {
    return area_leaf(height);
  }

  /* double competition_effect_state(Internals& vars) const { */
    /* return area_leaf_state(vars); */
  /* } */

  double compute_competition(double z, const Internals& vars) const {
    return area_leaf_above(z, vars.state(HEIGHT_INDEX), vars.state(AREA_LEAF_INDEX));
  }

  void refresh_indices();


  // ES20 Methods  ----------------------------------------------

  // [eqn 2] area_leaf (inverse of [eqn 3])
  double area_leaf(double height) const;

  // [eqn 1] mass_leaf (inverse of [eqn 2])
  double mass_leaf(double area_leaf) const;

  // [eqn 4] area and mass of sapwood
  double area_sapwood(double area_leaf) const;
  double mass_sapwood(double area_sapwood, double height) const;

  // [eqn 5] area and mass of bark
  double area_bark(double area_leaf) const;
  double mass_bark (double area_bark, double height) const;

  double area_stem(double area_bark, double area_sapwood,
                            double area_heartwood) const;
  double diameter_stem(double area_stem) const;

  // [eqn 7] Mass of (fine) roots
  double mass_root(double area_leaf) const;

  // [eqn 8] Total Mass
  double mass_live(double mass_leaf, double mass_bark,
                   double mass_sapwood, double mass_root) const;

  double mass_total(double mass_leaf, double mass_bark, double mass_sapwood,
                    double mass_heartwood, double mass_root) const;

  double mass_above_ground(double mass_leaf, double mass_bark,
                           double mass_sapwood, double mass_root) const;

  double dbiomass_dt(const ES20_Environment& environment, double mass_storage) const;
  
  double dmass_storage_dt(double net_mass_production_dt_, double dbiomass_dt_ ) const;
  
  double storage_portion(double mass_storage, double mass_leaf, double mass_bark, double mass_sapwood,
                         double mass_root) const;
  
  void compute_rates(const ES20_Environment& environment, bool reuse_intervals,
                Internals& vars);

  void update_dependent_aux(const int index, Internals& vars);


  // [eqn 13] Total maintenance respiration
  double respiration(double mass_leaf, double mass_sapwood,
                     double mass_bark, double mass_root, const ES20_Environment& environment) const;
  
  double respiration_leaf(double mass, const ES20_Environment& environment) const;
  double respiration_bark(double mass, const ES20_Environment& environment) const;
  double respiration_sapwood(double mass, const ES20_Environment& environment) const;
  double respiration_root(double mass, const ES20_Environment& environment) const;

  // [eqn 14] Total turnover
  double turnover(double mass_leaf, double mass_bark,
                  double mass_sapwood, double mass_root) const;
  double turnover_leaf(double mass) const;
  double turnover_bark(double mass) const;
  double turnover_sapwood(double mass) const;
  double turnover_root(double mass) const;

  // [eqn 15] Net production
  double net_mass_production_dt_A(double assimilation, double respiration) const;
  double net_mass_production_dt(const ES20_Environment& environment,
                                double height, double area_leaf, double mass_leaf, double mass_sapwood,
                                double mass_bark, double mass_root,
                                bool reuse_intervals=false);
  
  double net_mass_production_dt(const ES20_Environment& environment,
                                const Internals& vars, 
                                bool reuse_intervals=false);

  // [eqn 17] Rate of offspring production
  double fecundity_dt(double net_mass_production_dt,
                      double fraction_allocation_reproduction) const;

  // [eqn 18] Fraction of mass growth that is leaves
  double darea_leaf_dmass_live(double area_leaf, double height, double mass_sapwood) const;
  
  double mass_leaf_dt(double area_leaf, double darea_leaf_dt_) const; 
  // change in height per change in leaf area
  double dheight_darea_leaf(double area_leaf_, double height_) const;
  // adjustments made to allometry because of turnover time
  double height_area_leaf_allometry_difference(double area_leaf_, double height_) const;
  double height_area_leaf_allometry_adjustment(double area_leaf_, double height_) const;
  // Mass of leaf needed for new unit area leaf, d m_s / d a_l
  double dmass_leaf_darea_leaf(double area_leaf) const;
  // Mass of stem needed for new unit area leaf, d m_s / d a_l
  double dmass_sapwood_darea_leaf(double area_leaf, double height, double mass_sapwood) const;
  // adjustments made to allometry because of turnover time
  double mass_sapwood_area_leaf_allometry_difference(double area_leaf_, double height_, double mass_sapwood_) const;
  double mass_sapwood_area_leaf_allometry_adjustment(double area_leaf_, double height_, double mass_sapwood_) const;
  // Mass of bark needed for new unit area leaf, d m_b / d a_l
  double dmass_bark_darea_leaf(double area_leaf, double height, double mass_sapwood) const;
  // Mass of root needed for new unit area leaf, d m_r / d a_l
  double dmass_root_darea_leaf() const;
  // Growth rate of basal diameter_stem per unit stem area
  double ddiameter_stem_darea_stem(double area_stem) const;
  // Growth rate of components per unit time:
  double area_leaf_dt(double darea_leaf_dmass_live_, double dbiomass_dt_) const;
  double area_sapwood_dt(double area_leaf_dt) const;
  double mass_sapwood_dt(double dleaf_area_dt_, double dmass_sapwood_darea_leaf_) const;
  double area_heartwood_dt(double area_sapwood) const;
  double area_bark_dt(double area_leaf_dt) const;
  double mass_bark_dt(double dsapwood_mass_dt_) const;
  double area_stem_dt(double area_leaf, double area_leaf_dt) const;
  double diameter_stem_dt(double area_stem, double area_stem_dt) const;
  double mass_root_dt(double area_leaf_dt) const;
  double mass_live_dt(double fraction_allocation_reproduction,
                      double net_mass_production_dt) const;
  double mass_total_dt(double fraction_allocation_reproduction,
                       double net_mass_production_dt,
                       double mass_heartwood_dt) const;
  double mass_above_ground_dt(double area_leaf,
                              double fraction_allocation_reproduction,
                              double net_mass_production_dt,
                              double mass_heartwood_dt,
                              double area_leaf_dt) const;
  
  double mass_heartwood_dt(double mass_sapwood) const;
  
  double mass_live_given_height(double height) const;
  double height_given_mass_leaf(double mass_leaf_) const;
  
  
  double mortality_dt(double storage_portion, double productivity_area, const ES20_Environment& environment, double cumulative_mortality) const;
  double mortality_growth_independent_dt()const ;
  double mortality_storage_dependent_dt(double storage_portion) const;
  double mortality_growth_dependent_dt(double productivity_area, const ES20_Environment& environment) const;
  // [eqn 20] Survival of seedlings during establishment
  double establishment_probability(const ES20_Environment& environment);
  
  // * Competitive environment
  // [eqn 11] total leaf area above height above height `z` for given plant
  double area_leaf_above(double z, double height, double area_leaf) const;
  // [eqn 10] Fraction of leaf area above height `z`
  double Q(double z, double height) const;

  // The aim is to find a plant height that gives the correct seed mass.
  double height_seed(void) const;

  // Set constants within ES20_Strategy
  void prepare_strategy();

  // Initialize all relevant states
  void initialize_states(Internals& vars);

  // * Core traits
  double lma, rho, hmat, omega;
  // * Individual allometry
  // Canopy shape parameters
  double eta, eta_c;
  // Sapwood area per leaf area
  double theta;
  // Empirical constants for scaling relationships
  double a_l1, a_l2, a_r1;
  // Bark area per sapwood area
  double a_b1;
  // * Production
  // Respiration constants
  double r_s, r_b, r_r, r_l;
  // Yield = carbon fixed in tissue per carbon assimilated;
  double a_y;
  // Conversion factor
  double a_bio;
  // Leaf, bark sapwood, and root turnover rates
  double k_l, k_b, k_s, k_r;
  // Leaf productivity parameters  - only used when no N reallocation
  double a_p1, a_p2;

  // * Seed production
  // Accessory cost of reproduction, per seed
  double a_f3;
  // Proportion production allocated to reproduction
  double a_f1;
  // Size range across which individuals mature
  double a_f2;

  // * Mortality
  // Probability of survival during dispersal
  double S_D;
  // Parameter for seedling mortality
  double a_d0;
  // Baseline structural mortality rate
  double d_I;
  // Baseline for storage mortality rate
  double a_dS1;
  // Coefficient for storage production in mortality function
  double a_dS2;
 // Baseline for growth mortality rate
  double a_dG1;
  // Coefficient for dry mass production in mortality function
  double a_dG2;

  // * Storage
  // Time of growth switch (day)
  double t_s;
  // Proportion of storage allocated to biomass
  double a_s;
  // Proportion of storage as related to initial biomass
  double b_s1;
  // Speed of growth switch
  double b_s2;
  
  // Height and leaf area of a (germinated) seed
  double height_0;
  double area_leaf_0;
  double mass_leaf_0;
  double mass_sapwood_0;
  double mass_bark_0;
  double mass_root_0;
  double area_sapwood_0;
  double area_bark_0;
  double area_stem_0;
  double mass_storage_0;
  double diameter_stem_0;

  std::string name;

  ES20_Assimilation assimilator;

  
};

ES20_Strategy::ptr make_strategy_ptr(ES20_Strategy s);

}

#endif
