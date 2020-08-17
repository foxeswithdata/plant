#include <memory>
#include <plant/control.h>
#include <plant/qag_internals.h> // quadrature::intervals_type
#include <plant/internals.h> // quadrature::intervals_type
#include <plant/strategy.h>
#include <plant/models/es20_environment.h>
#include <plant/models/es20_strategy.h>
#include <plant/models/assimilation.h>
#include <plant/uniroot.h>

#include <RcppCommon.h> // NA_REAL

#include <iostream>

namespace plant {

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?
ES20_Strategy::ES20_Strategy() {
  // * Core traits - default values
  lma       = 0.1978791;  // Leaf mass per area [kg / m2]
  rho       = 608.0;      // Wood density [kg/m3]
  hmat      = 16.5958691; // Height at maturation [m]
  omega     = 3.8e-5;     // Seed mass [kg]
  
  // * Individual allometry
  // Canopy shape parameter (extra calculation here later)
  eta       = 12.0; // [dimensionless]
  // Ratio sapwood area area to leaf area
  theta     = 1.0/4669; // [dimensionless]
  // Height - leaf mass scaling
  // a_l1        = 5.44; // height with 1m2 leaf [m]
  // a_l2        = 0.306; // dimensionless scaling of height with leaf area
  // Euc Saligna
  a_l1 = 1.585053;
  a_l2 = 0.4804506;
  
  // Root mass per leaf area
  a_r1        = 0.07;  //[kg / m]
  // Ratio of bark area : sapwood area
  a_b1         = 0.17; // [dimensionless]
  
  // * Production
  // Ratio of leaf dark respiration to leaf mass [mol CO2 / yr  / kg]
  // =  [mol CO2 / m2 / yr]  |  (39.27 = 2100 * 0.00187)  | narea * photosynthesis_per_nitrogen
  //    / [kg(leaf) / m2 ]   |    / (0.1978791)           | lma
  // Hard coded in value of lma here so that this value doesn't change
  // if that trait changes above.
  r_l   = 39.27 / 0.1978791;
  // Root respiration per mass [mol CO2 / yr / kg]
  r_r   = 217.0;
  // Sapwood respiration per stem mass  [mol CO2 / yr / kg]
  // = respiration per volume [mol CO2 / m3 / yr]
  // /  wood density [kg/m3]
  r_s   = 4012.0 / 608.0;
  // Bark respiration per stem mass
  // assumed to be twice rate of sapwood
  // (NOTE that there is a re-parametrisation here relative to the paper
  // -- r_b is defined (new) as 2*r_s, whereas the paper assumes a
  // fixed multiplication by 2)
  r_b   = 2.0 * r_s;
  // Carbon conversion parameter
  a_y      = 0.7;
  // Constant converting assimilated CO2 to dry mass [kg / mol]
  // (12E-3 / 0.49)
  a_bio  = 2.45e-2;
  // Leaf turnover [/yr]
  k_l=  0.4565855;
  // Bark turnover [/yr]
  k_b    = 0.2;
  // Sapwood turnover [/yr]
  k_s     = 0.2;
  // Root turnover [/yr]
  k_r    = 1.0;
  // Parameters of the hyperbola for annual LRC
  a_p1   = 151.177775377968; // [mol CO2 / yr / m2]
  a_p2   = 0.204716166503633; // [dimensionless]
  
  // * Seed production
  // Accessory cost of reproduction
  a_f3  = 3.0 *  3.8e-5; // [kg per seed]
  
  // Maximum allocation to reproduction
  a_f1   = 0.05; //[dimensionless]
  // Size range across which individuals mature
  a_f2   = 50; // [dimensionless]
  
  // * Mortality parameters
  // Probability of survival during dispersal
  S_D   = 0.25; // [dimensionless]
  // Parameter for seedling survival
  a_d0    = 0.1; //[kg / yr / m2]
  // Baseline for intrinsic mortality
  d_I    = 0.01; // [ / yr]
  // Baseline rate for growth-related mortality
  a_dG1    = 5.5; // [ / yr]
  // Risk coefficient for dry mass production (per area)
  a_dG2    = 20.0;// [yr m2 / kg ]
  
  
  // storage parameters
  // storage allocation parameter [kg kg-1]
  // a_s = 0.06;
  a_s = 0.06 * 365;
  // time of switch [yr]
  t_s = 0.4109589;
  
  
  // Will get computed properly by prepare_strategy
  height_0 = NA_REAL;
  area_leaf_0 = NA_REAL;
  mass_leaf_0 = NA_REAL;
  mass_sapwood_0 = NA_REAL;
  mass_bark_0 = NA_REAL;
  mass_root_0 = NA_REAL;
  area_sapwood_0 = NA_REAL;
  area_bark_0 = NA_REAL;
  area_stem_0 = NA_REAL;
  mass_storage_0 = NA_REAL;
  eta_c    = NA_REAL;
  
  collect_all_auxillary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "ES20";
}

void ES20_Strategy::refresh_indices () {
  // Create and fill the name to state index maps
  state_index = std::map<std::string,int>();
  aux_index   = std::map<std::string,int>();
  std::vector<std::string> aux_names_vec = aux_names();
  std::vector<std::string> state_names_vec = state_names();
  for (int i = 0; i < state_names_vec.size(); i++) {
    state_index[state_names_vec[i]] = i;
  }
  for (int i = 0; i < aux_names_vec.size(); i++) {
    aux_index[aux_names_vec[i]] = i;
  }
}


// [eqn 1] mass_leaf (inverse of [eqn 2])
double ES20_Strategy::mass_leaf(double area_leaf) const {
  return area_leaf * lma;
}

// [eqn 4] area and mass of sapwood
double ES20_Strategy::area_sapwood(double area_leaf) const {
  return area_leaf * theta;
}

double ES20_Strategy::mass_sapwood(double area_sapwood, double height) const {
  return area_sapwood * height * eta_c * rho;
}

// [eqn 5] area and mass of bark
double ES20_Strategy::area_bark(double area_leaf) const {
  return a_b1 * area_leaf * theta;
}

double ES20_Strategy::mass_bark(double area_bark, double height) const {
  return area_bark * height * eta_c * rho;
}

double ES20_Strategy::area_stem(double area_bark, double area_sapwood,
                                double area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

double ES20_Strategy::diameter_stem(double area_stem) const {
  return std::sqrt(4 * area_stem / M_PI);
}

// [eqn 7] Mass of (fine) roots
double ES20_Strategy::mass_root(double area_leaf) const {
  return a_r1 * area_leaf;
}

// [eqn 8] Total mass
double ES20_Strategy::mass_live(double mass_leaf, double mass_bark,
                                double mass_sapwood, double mass_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root;
}

double ES20_Strategy::mass_total(double mass_leaf, double mass_bark,
                                 double mass_sapwood, double mass_heartwood,
                                 double mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root;
}

double ES20_Strategy::mass_above_ground(double mass_leaf, double mass_bark,
                                        double mass_sapwood, double mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_root;
}

// for updating auxillary state
void ES20_Strategy::update_dependent_aux(const int index, Internals& vars) {
  if (index == HEIGHT_INDEX) {
    double height = vars.state(HEIGHT_INDEX);
    vars.set_aux(aux_index.at("competition_effect"), area_leaf(height));
  }
}

double ES20_Strategy::dbiomass_dt(const ES20_Environment& environment, 
                                  double mass_storage) const {
  if (environment.time_in_year() < t_s){
    return mass_storage * a_s;
  }
  else {
    return 0;
  }
}

double ES20_Strategy::dmass_storage_dt(double net_mass_production_dt_, double dbiomass_dt_) const {
  return net_mass_production_dt_ - dbiomass_dt_;
}


// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
void ES20_Strategy::compute_rates(const ES20_Environment& environment,
                                  bool reuse_intervals,
                                  Internals& vars) {
  
  double height = vars.state(HEIGHT_INDEX);
  double area_leaf_ = vars.state(AREA_LEAF_INDEX);
  double area_bark_ = vars.state(AREA_BARK_INDEX);
  double area_heartwood_ = vars.state(AREA_HEARTWOOD_INDEX);
  double area_sapwood_ = vars.state(AREA_SAPWOOD_INDEX);
  double mass_leaf_ = vars.state(MASS_LEAF_INDEX);
  double mass_bark_ = vars.state(MASS_BARK_INDEX);
  double mass_sapwood_ = vars.state(MASS_SAPWOOD_INDEX);
  double mass_root_ = vars.state(MASS_ROOT_INDEX);
  double mass_storage_ = vars.state(MASS_STORAGE_INDEX);
  
  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height, area_leaf_, mass_leaf_, mass_sapwood_,
                           mass_bark_, mass_root_, reuse_intervals);
  
  const double dbiomass_dt_ = dbiomass_dt(environment, mass_storage_);
  double dmass_storage_dt_  = dmass_storage_dt(net_mass_production_dt_, dbiomass_dt_);
  
  // store the aux sate
  vars.set_aux(aux_index.at("net_mass_production_dt"), net_mass_production_dt_);
  vars.set_aux(aux_index.at("dbiomass_dt"), dbiomass_dt_);
  vars.set_aux(aux_index.at("respiration_leaf_dt"), a_bio * a_y * respiration_leaf(mass_leaf_,environment));
  vars.set_aux(aux_index.at("respiration_bark_dt"), a_bio * a_y * respiration_bark(mass_bark_,environment));
  vars.set_aux(aux_index.at("respiration_root_dt"), a_bio * a_y * respiration_root(mass_root_,environment));
  vars.set_aux(aux_index.at("respiration_sapwood_dt"), a_bio * a_y * respiration_sapwood(mass_sapwood_,environment));
  vars.set_aux(aux_index.at("respiration_dt"), a_bio * a_y * (respiration_leaf(mass_leaf_,environment) +
    respiration_bark(mass_bark_,environment) + respiration_sapwood(mass_sapwood_,environment) + respiration_root(mass_root_,environment)));
  
  if (dbiomass_dt_ > 0) {
    
    // Changes in height and leaf area 
    
    const double dheight_darea_leaf_ = dheight_darea_leaf(area_leaf_);
    
    const double darea_leaf_dmass_live_ = darea_leaf_dmass_live(area_leaf_, height);
    const double darea_leaf_dt_ = darea_leaf_dmass_live_ * dbiomass_dt_;
    const double dmass_leaf_dt_ = mass_leaf_dt(area_leaf_, darea_leaf_dt_);
    
    vars.set_aux(aux_index.at("area_leaf_a_l_dt"), area_leaf_ + darea_leaf_dt_);
    
    const double dheight_dt_ = dheight_darea_leaf_ * darea_leaf_dt_;
    
    vars.set_aux(aux_index.at("darea_leaf_dmass_live"), darea_leaf_dmass_live_);
    
    
    vars.set_rate(HEIGHT_INDEX, dheight_dt_);
    vars.set_rate(AREA_LEAF_INDEX, darea_leaf_dt_);
    vars.set_rate(MASS_LEAF_INDEX, dmass_leaf_dt_ - turnover_leaf(mass_leaf_));
    
    // Changes in sapwood and heartwood
    
    const double darea_heartwood_dt_ = area_heartwood_dt(area_sapwood_);
    const double dmass_heartwood_dt_ = mass_heartwood_dt(mass_sapwood_);
    
    vars.set_rate(state_index.at("area_heartwood"), darea_heartwood_dt_);
    vars.set_rate(state_index.at("mass_heartwood"), dmass_heartwood_dt_);
    
    const double darea_sapwood_dt_ = area_sapwood_dt(darea_leaf_dt_);
    const double dmass_sapwood_darea_leaf_ = dmass_sapwood_darea_leaf(area_leaf_, height);
    const double dmass_sapwood_dt_ = mass_sapwood_dt(darea_leaf_dt_, dmass_sapwood_darea_leaf_);
    
    vars.set_rate(state_index.at("area_sapwood"), darea_sapwood_dt_ - turnover_sapwood(area_sapwood_));
    vars.set_rate(state_index.at("mass_sapwood"), dmass_sapwood_dt_ - turnover_sapwood(mass_sapwood_));
    
    // changes in bark
    const double darea_bark_dt_ = area_bark_dt(darea_leaf_dt_);
    const double dmass_bark_dt_ = mass_bark_dt(dmass_sapwood_dt_);
    
    vars.set_rate(state_index.at("area_bark"), darea_bark_dt_ - turnover_bark(area_bark_));
    vars.set_rate(state_index.at("mass_bark"), dmass_bark_dt_ - turnover_bark(mass_bark_));
    
    // changes in stem
    const double darea_stem_dt_ = darea_sapwood_dt_ + darea_heartwood_dt_ + darea_bark_dt_ - turnover_bark(area_bark_);
    const double ddiameter_stem_dt_ = diameter_stem_dt(area_stem(area_heartwood_, area_sapwood_, area_bark_), darea_stem_dt_);
    
    vars.set_rate(state_index.at("diameter_stem"), ddiameter_stem_dt_);
    vars.set_rate(state_index.at("area_stem"), darea_stem_dt_);
    
    // changes in root
    const double dmass_root_dt_ = mass_root_dt(darea_leaf_dt_);
    
    vars.set_rate(state_index.at("mass_root"), dmass_root_dt_ - turnover_root(mass_root_));
    
    // fecundity 
    vars.set_rate(FECUNDITY_INDEX,
                  fecundity_dt(mass_storage_, height));
    
    // update losses due to reproduction
    dmass_storage_dt_ = dmass_storage_dt_ - fecundity_dt(mass_storage_, height);
    
  } else {
    vars.set_rate(HEIGHT_INDEX, 0.0);
    vars.set_rate(FECUNDITY_INDEX, 0.0);
    vars.set_rate(AREA_LEAF_INDEX, - turnover_leaf(area_leaf_));
    vars.set_rate(MASS_LEAF_INDEX, - turnover_leaf(mass_leaf_));
    vars.set_rate(state_index.at("area_heartwood"), turnover_sapwood(area_sapwood_));
    vars.set_rate(state_index.at("mass_heartwood"), turnover_sapwood(mass_sapwood_));
    vars.set_rate(state_index.at("area_sapwood"), - turnover_sapwood(area_sapwood_));
    vars.set_rate(state_index.at("mass_sapwood"), (- turnover_sapwood(mass_sapwood_))); // <- error in this (?)
    vars.set_rate(state_index.at("area_bark"), - turnover_bark(area_bark_));
    vars.set_rate(state_index.at("mass_bark"), - turnover_bark(mass_bark_));
    vars.set_rate(state_index.at("diameter_stem"), 0.0);
    vars.set_rate(state_index.at("mass_root"), - turnover_root(mass_root_));
    vars.set_rate(state_index.at("area_stem"), - turnover_bark(area_bark_));
  }
  // [eqn 21] - Instantaneous mortality rate
  vars.set_rate(MASS_STORAGE_INDEX, dmass_storage_dt_);
  vars.set_rate(MORTALITY_INDEX,
                mortality_dt(mass_storage_ / mass_live(mass_leaf_, mass_bark_, mass_sapwood_, mass_root_), vars.state(MORTALITY_INDEX)));
}


// [eqn 13] Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by a_y*a_bio.
double ES20_Strategy::respiration(double mass_leaf, double mass_sapwood,
                                  double mass_bark, double mass_root, const ES20_Environment& environment) const {
  return respiration_leaf(mass_leaf, environment) +
    respiration_bark(mass_bark,environment) +
    respiration_sapwood(mass_sapwood,environment) +
    respiration_root(mass_root, environment);
}

double ES20_Strategy::respiration_leaf(double mass, const ES20_Environment& environment) const {
  if (environment.stressed()){
    return (r_l/2) * mass;
  }
  return r_l * mass;
}

double ES20_Strategy::respiration_bark(double mass, const ES20_Environment& environment) const {
  if (environment.stressed()){
    return (r_b/2) * mass;
  }
  return r_b * mass;
}

double ES20_Strategy::respiration_sapwood(double mass, const ES20_Environment& environment) const {
  if (environment.stressed()){
    return (r_s/2) * mass;
  }
  return r_s * mass;
}

double ES20_Strategy::respiration_root(double mass, const ES20_Environment& environment) const {
  if (environment.stressed()){
    return (r_r/2) * mass;
  }
  return r_r * mass;
}

// [eqn 14] Total turnover
double ES20_Strategy::turnover(double mass_leaf, double mass_bark,
                               double mass_sapwood, double mass_root) const {
  return turnover_leaf(mass_leaf) +
    turnover_bark(mass_bark) +
    turnover_sapwood(mass_sapwood) +
    turnover_root(mass_root);
}

double ES20_Strategy::turnover_leaf(double mass) const {
  return 0;
  // return k_l * mass;
}

double ES20_Strategy::turnover_bark(double mass) const {
  return 0;
  // return k_b * mass;
}

double ES20_Strategy::turnover_sapwood(double mass) const {
  return 0;
  // return k_s * mass;
}

double ES20_Strategy::turnover_root(double mass) const {
  return 0;
  // return k_r * mass;
}

// [eqn 15] Net production
//
// NOTE: Translation of variable names from the Falster 2011.  Everything
// before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
double ES20_Strategy::net_mass_production_dt_A(double assimilation, double respiration) const {
  return a_bio * a_y * (assimilation - respiration);
}

// One shot calculation of net_mass_production_dt
// Used by establishment_probability() and compute_rates().
double ES20_Strategy::net_mass_production_dt(const ES20_Environment& environment,
                                             double height, double area_leaf, double mass_leaf, double mass_sapwood,
                                             double mass_bark, double mass_root,
                                             bool reuse_intervals) {
  const double assimilation_ = assimilator.assimilate(control, environment, height,
                                                      area_leaf, reuse_intervals);
  const double respiration_ =
    respiration(mass_leaf, mass_sapwood, mass_bark, mass_root, environment);
  return net_mass_production_dt_A(assimilation_, respiration_);
}


double ES20_Strategy::area_leaf_dt(double darea_leaf_dmass_live_, double dbiomass_dt_) const{
  return darea_leaf_dmass_live_ * dbiomass_dt_;
}

double ES20_Strategy::mass_leaf_dt(double area_leaf, double darea_leaf_dt_) const {
  return dmass_leaf_darea_leaf(area_leaf) * darea_leaf_dt_ ;
}

// [eqn 2] area_leaf (inverse of [eqn 3])
double ES20_Strategy::area_leaf(double height) const {
  return pow(height / a_l1, 1.0 / a_l2);
}

// [eqn 17] Rate of offspring production
double ES20_Strategy::fecundity_dt(double mass_storage, double height) const {
  return mass_storage * a_f1 / ((omega + a_f3) * (1.0 + exp(a_f2 * (1.0 - height / hmat))));
}

double ES20_Strategy::darea_leaf_dmass_live(double area_leaf, double height) const {
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
                  + dmass_sapwood_darea_leaf(area_leaf, height)
                  + dmass_bark_darea_leaf(area_leaf, height)
                  + dmass_root_darea_leaf());
}

// TODO: Ordering below here needs working on, probably as @dfalster
// does equation documentation?
double ES20_Strategy::dheight_darea_leaf(double area_leaf) const {
  return a_l1 * a_l2 * pow(area_leaf, a_l2 - 1);
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
double ES20_Strategy::dmass_leaf_darea_leaf(double /* area_leaf */) const {
  return lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
double ES20_Strategy::dmass_sapwood_darea_leaf(double area_leaf, double height) const {
  return rho * eta_c *theta * (height + a_l1 * a_l2 * pow(area_leaf, a_l2));
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
double ES20_Strategy::dmass_bark_darea_leaf(double area_leaf, double height) const {
  return a_b1 * dmass_sapwood_darea_leaf(area_leaf, height);
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
double ES20_Strategy::dmass_root_darea_leaf() const {
  return a_r1;
}

// Growth rate of basal diameter_stem per unit time
double ES20_Strategy::ddiameter_stem_darea_stem(double area_stem) const {
  return pow(M_PI * area_stem, -0.5);
}

// Growth rate of sapwood area at base per unit time
double ES20_Strategy::area_sapwood_dt(double area_leaf_dt) const {
  return area_leaf_dt * theta;
}

double ES20_Strategy::mass_sapwood_dt(double dleaf_area_dt_, double dmass_sapwood_darea_leaf_) const{
  return dleaf_area_dt_ * dmass_sapwood_darea_leaf_;
}

// Note, unlike others, heartwood growth does not depend on leaf area growth, but
// rather existing sapwood
double ES20_Strategy::area_heartwood_dt(double area_sapwood) const {
  return k_s * area_sapwood;
}

// Growth rate of bark area at base per unit time
double ES20_Strategy::area_bark_dt(double area_leaf_dt) const {
  return a_b1 * area_leaf_dt * theta;
}

// Growth rate of bark mass at base per unit time
double ES20_Strategy::mass_bark_dt(double dsapwood_mass_dt) const {
  return a_b1 * dsapwood_mass_dt;
}

// Growth rate of stem basal area per unit time
double ES20_Strategy::area_stem_dt(double area_leaf,
                                    double area_leaf_dt) const {
  return area_sapwood_dt(area_leaf_dt) +
    area_bark_dt(area_leaf_dt) +
    area_heartwood_dt(area_leaf);
}

// Growth rate of basal diameter_stem per unit time
double ES20_Strategy::diameter_stem_dt(double area_stem, double area_stem_dt) const {
  return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
}

// Growth rate of root mass per unit time
double ES20_Strategy::mass_root_dt(double area_leaf_dt) const {
  return area_leaf_dt * dmass_root_darea_leaf();
}

double ES20_Strategy::mass_live_dt(double fraction_allocation_reproduction,
                                   double net_mass_production_dt) const {
  return (1 - fraction_allocation_reproduction) * net_mass_production_dt;
}

// TODO: Change top two to use mass_live_dt
double ES20_Strategy::mass_total_dt(double fraction_allocation_reproduction,
                                    double net_mass_production_dt,
                                    double mass_heartwood_dt) const {
  return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
    mass_heartwood_dt;
}

// TODO: Do we not track root mass change?
double ES20_Strategy::mass_above_ground_dt(double area_leaf,
                                           double fraction_allocation_reproduction,
                                           double net_mass_production_dt,
                                           double mass_heartwood_dt,
                                           double area_leaf_dt) const {
  const double mass_root_dt =
    area_leaf_dt * dmass_root_darea_leaf();
  return mass_total_dt(fraction_allocation_reproduction, net_mass_production_dt,
                       mass_heartwood_dt) - mass_root_dt;
}

double ES20_Strategy::mass_heartwood_dt(double mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


double ES20_Strategy::height_given_mass_leaf(double mass_leaf) const {
  return a_l1 * pow(mass_leaf / lma, a_l2);
}

double ES20_Strategy::mortality_dt(double storage_portion,
                                   double cumulative_mortality) const {
  
  // NOTE: When plants are extremely inviable, the rate of change in
  // mortality can be Inf, because net production is negative, leaf
  // area is small and so we get exp(big number).  However, most of
  // the time that happens we should get infinite mortality variable
  // levels and the rate of change won't matter.  It is possible that
  // we will need to trim this to some large finite value, but for
  // now, just checking that the actual mortality rate is finite.
  if (R_FINITE(cumulative_mortality)) {
    return
    mortality_growth_independent_dt() +
      mortality_growth_dependent_dt(storage_portion);
  } else {
    // If mortality probability is 1 (latency = Inf) then the rate
    // calculations break.  Setting them to zero gives the correct
    // behaviour.
    return 0.0;
  }
}

double ES20_Strategy::mortality_growth_independent_dt() const {
  return d_I;
}

double ES20_Strategy::mortality_growth_dependent_dt(double storage_portion) const {
  if(storage_portion <= 0){
    return 1;
  }
  return a_dG1 * exp(-a_dG2 * storage_portion);
}

// [eqn 20] Survival of seedlings during establishment
double ES20_Strategy::establishment_probability(const ES20_Environment& environment) {
  
  //TODO: add dependency on stress
  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height_0, area_leaf_0,
                           mass_leaf_0, mass_sapwood_0,
                           mass_bark_0, mass_root_0);
  if (net_mass_production_dt_ > 0) {
    const double tmp = a_d0 * area_leaf_0 / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0);
  } else {
    return 0.0;
  }
}

double ES20_Strategy::area_leaf_above(double z, double height) const {
  return area_leaf(height) * Q(z, height);
}

// [eqn 10] ... Fraction of leaf area above height 'z' for an
//              individual of height 'height'
double ES20_Strategy::Q(double z, double height) const {
  if (z > height) {
    return 0.0;
  }
  const double tmp = 1.0-pow(z / height, eta);
  return tmp * tmp;
}

double ES20_Strategy::mass_live_given_height(double height) const {
  double area_leaf_ = area_leaf(height);
  return mass_leaf(area_leaf_) +
    mass_bark(area_bark(area_leaf_), height) +
    mass_sapwood(area_sapwood(area_leaf_), height) +
    mass_root(area_leaf_);
}

// The aim is to find a plant height that gives the correct seed mass.
double ES20_Strategy::height_seed(void) const {
  
  // Note, these are not entirely correct bounds. Ideally we would use height
  // given *total* mass, not leaf mass, but that is difficult to calculate.
  // Using "height given leaf mass" will expand upper bound, but that's ok
  // most of time. Only issue is that could break with obscure parameter
  // values for LMA or height-leaf area scaling. Could instead use some
  // absolute maximum height for new seedling, e.g. 1m?
  const double
  h0 = height_given_mass_leaf(std::numeric_limits<double>::min()),
    h1 = height_given_mass_leaf(omega);
  
  const double tol = control.plant_seed_tol;
  const size_t max_iterations = control.plant_seed_iterations;
  
  auto target = [&] (double x) mutable -> double {
    return mass_live_given_height(x) - omega;
  };
  
  // return util::uniroot(target, h0, h1, tol, max_iterations);
  
  return 0.4;
}

void ES20_Strategy::prepare_strategy() {
  // Set up the integrator
  control.initialize();
  assimilator.initialize(a_p1, a_p2, eta);
  // NOTE: this pre-computes something to save a very small amount of time
  eta_c = 1 - 2/(1 + eta) + 1/(1 + 2*eta);
  // NOTE: Also pre-computing, though less trivial
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);
  mass_leaf_0 = mass_leaf(area_leaf_0);
  area_sapwood_0 = area_sapwood(area_leaf_0);
  mass_sapwood_0 = mass_sapwood(area_sapwood_0, height_0);
  area_bark_0 = area_bark(area_leaf_0);
  mass_bark_0 = mass_bark(area_bark_0, height_0);
  mass_root_0 = mass_root(area_leaf_0);
  area_stem_0 = area_stem(area_bark_0, area_sapwood_0, 0); 
  mass_storage_0 = 0.1 * mass_live(mass_leaf_0, mass_bark_0, mass_sapwood_0, mass_root_0);
  diameter_stem_0 = diameter_stem(area_stem_0);
}

void ES20_Strategy::initialize_states(Internals &vars){
  vars.set_state(HEIGHT_INDEX, height_0);
  vars.set_state(AREA_LEAF_INDEX, area_leaf_0);
  vars.set_state(MASS_LEAF_INDEX, mass_leaf_0);
  vars.set_state(state_index.at("area_sapwood"), area_sapwood_0);
  vars.set_state(state_index.at("mass_sapwood"), mass_sapwood_0);
  vars.set_state(state_index.at("area_bark"), area_bark_0);
  vars.set_state(state_index.at("mass_bark"), mass_bark_0);
  vars.set_state(state_index.at("mass_root"), mass_root_0);
  vars.set_state(state_index.at("area_stem"), area_stem_0);
  vars.set_state(state_index.at("mass_storage"), mass_storage_0);
  vars.set_state(state_index.at("diameter_stem"), diameter_stem_0);
}

ES20_Strategy::ptr make_strategy_ptr(ES20_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<ES20_Strategy>(s);
}
}
