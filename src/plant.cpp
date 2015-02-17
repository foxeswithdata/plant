#include <tree2/plant.h>
#include <Rcpp.h>
#include <functional>

namespace tree2 {

Plant::Plant(strategy_ptr_type s)
  : strategy(s) {
  set_height(strategy->height_0);
}

// Individual size
// [eqn 1-8] Update size variables given an input leaf mass

// Height is really important -- everything else follows from it.
double Plant::height() const {
  return vars.height;
}
double Plant::height_rate() const {
  return vars.height_growth_rate;
}
// NOTE: Only recomputes size variables if the height is actually
// different.  This is totally safe if nothing else sets either height
// or any size variable except this method.  This could help save
// quite a bit of calculation time and book-keeping down the track.
// If we're off because of a floating point difference, the worst that
// happens is that we recompute the variables again.
void Plant::set_height(double height_) {
  if (height_ < 0.0) {
    Rcpp::stop("height must be positive (given " +
               util::to_string(height_) + ")");
  }
  // TODO: in the original version of the model, all of plant size was driven
  // only  by height, so we only needed to check this here.  But now we have
  // two size variables (heartwood mass being the other) we need to be more
  // careful. For the new plant type, replace this height function with a new
  // function that takes both size variables as arguments and update all the
  // size variables at that point.

  if (!util::identical(height_, height())) {
    compute_vars_size(height_);
  }
}

double Plant::mortality() const {
  return vars.mortality;
}
double Plant::mortality_rate() const {
  return vars.mortality_rate;
}
void Plant::set_mortality(double x) {
  vars.mortality = x;
}

double Plant::fecundity() const {
  return vars.fecundity;
}
double Plant::fecundity_rate() const {
  return vars.fecundity_rate;
}
void Plant::set_fecundity(double x) {
  vars.fecundity = x;
}

double Plant::heartwood_area() const {
  return vars.heartwood_area;
}

double Plant::heartwood_area_rate() const {
  return vars.heartwood_area_rate;
}

void Plant::set_heartwood_area(double x) {
  // TODO: Consider recomputing the size variables here
  // See set_heartwood_mass
  vars.heartwood_area = x;
}

double Plant::heartwood_mass() const {
  return vars.heartwood_mass;
}

double Plant::heartwood_mass_rate() const {
  return vars.heartwood_mass_rate;
}

void Plant::set_heartwood_mass(double x) {
  // TODO: This needs to update size variables but does not yet,
  // and won't until we get the new version done.
  // See notes in set_height.
  vars.heartwood_mass = x;
}

// * Competitive environment
double Plant::leaf_area() const {
  return vars.leaf_area;
}

// [      ] Leaf area (not fraction) above height `z`
double Plant::leaf_area_above(double z) const {
  if (z < 0.0) {
    Rcpp::stop("Negative heights do not make sense");
  }
  return strategy->leaf_area_above(z, vars.height, vars.leaf_area);
}

// * Mass production
// [eqn 12-19,21] Update physiological variables given the current
// light environment (and given the current set of size variables).
void Plant::compute_vars_phys(const Environment& environment,
			      bool reuse_intervals) {
  // [eqn 12] Gross annual CO2 assimilation
  vars.assimilation = strategy->assimilation(environment, vars.height,
                                             vars.leaf_area, reuse_intervals);

  // [eqn 13] Total maintenance respiration
  vars.respiration = strategy->respiration(vars.leaf_mass,
                                           vars.sapwood_mass,
                                           vars.bark_mass,
                                           vars.root_mass);

  // [eqn 14] Total turnover
  vars.turnover = strategy->turnover(vars.leaf_mass, vars.sapwood_mass,
                                     vars.bark_mass, vars.root_mass);

  // [eqn 15] Net production
  vars.net_mass_production = strategy->net_mass_production(vars.assimilation,
                                                 vars.respiration,
                                                 vars.turnover);

  if (vars.net_mass_production > 0) {
    // [eqn 16] - Fraction of whole plant growth that is leaf
    vars.reproduction_mass_fraction =
      strategy->reproduction_mass_fraction(vars.height);

    // [eqn 17] - Rate of offspring production
    //
    // NOTE: In EBT, was multiplied by Pi_0 (survival during
    // dispersal), but we do not here.
    // NOTE: This is also a hyperparametrisation and should move into
    // the initialisation function.
    vars.fecundity_rate =
      strategy->dfecundity_dt(vars.net_mass_production,
                              vars.reproduction_mass_fraction);

    // [eqn 19] - Growth rate in leaf height
    // different to Falster 2010, which was growth rate in leaf mass
    vars.leaf_area_deployment_mass =
      strategy->leaf_area_deployment_mass(vars.leaf_area);
    vars.growth_mass_fraction = strategy->growth_mass_fraction(vars.height);

    vars.leaf_area_growth_rate =
      vars.net_mass_production * vars.growth_mass_fraction *
      vars.leaf_area_deployment_mass;
   vars.height_growth_rate =
      strategy->dheight_dleaf_area(vars.leaf_area) *
      vars.leaf_area_growth_rate;
   vars.heartwood_area_rate =
     strategy->dheartwood_area_dt(vars.leaf_area);
   vars.heartwood_mass_rate =
      strategy->dheartwood_mass_dt(vars.sapwood_mass);
  } else {
    vars.reproduction_mass_fraction = 0.0;
    vars.growth_mass_fraction       = 0.0;
    vars.fecundity_rate        = 0.0;
    vars.leaf_area_growth_rate = 0.0;
    vars.leaf_area_deployment_mass = 0.0;
    vars.height_growth_rate    = 0.0;
    vars.heartwood_area_rate   = 0.0;
    vars.heartwood_mass_rate   = 0.0;
  }

  // [eqn 21] - Instantaneous mortality rate
  //
  // NOTE: When plants are extremely inviable, the rate of change in
  // mortality can be Inf, because net production is negative, leaf
  // area is small and so we get exp(big number).  However, most of
  // the time that happens we should get infinite mortality variable
  // levels and the rate of change won't matter.  It is possible that
  // we will need to trim this to some large finite value, but for
  // now, just checking that the actual mortality rate is finite.
  if (R_FINITE(vars.mortality)) {
    vars.mortality_rate =
      strategy->mortality_dt(vars.net_mass_production / vars.leaf_area);
  } else {
    // If mortality probability is 1 (latency = Inf) then the rate
    // calculations break.  Setting them to zero gives the correct
    // behaviour.
    vars.mortality_rate = 0.0;
  }
}

// Extra accounting.
// TODO: This will move into the "super size" plant.
void Plant::compute_vars_growth() {
  const strategy_type *s = strategy.get(); // for brevity.
  const double leaf_area = vars.leaf_area,
    leaf_area_growth_rate = vars.leaf_area_growth_rate;

  // Changes with leaf area:
  vars.dheight_dleaf_area       = s->dheight_dleaf_area(leaf_area);
  vars.dsapwood_mass_dleaf_area = s->dsapwood_mass_dleaf_area(leaf_area);
  vars.dbark_mass_dleaf_area    = s->dbark_mass_dleaf_area(leaf_area);
  vars.droot_mass_dleaf_area    = s->droot_mass_dleaf_area(leaf_area);

  // Changes over time:
  vars.dsapwood_area_dt    = s->dsapwood_area_dt(leaf_area_growth_rate);
  vars.dbark_area_dt       = s->dbark_area_dt(leaf_area_growth_rate);
  vars.dstem_area_dt      = s->dstem_area_dt(leaf_area,
                                               leaf_area_growth_rate);
  vars.dstem_diameter_dt      = s->dstem_diameter_dt(leaf_area,
                                          vars.bark_area,
                                          vars.sapwood_area,
                                          vars.heartwood_area,
                                          leaf_area_growth_rate);
  vars.droot_mass_dt       = s->droot_mass_dt(leaf_area,
                                         leaf_area_growth_rate);
  vars.dlive_mass_dt       = s->dlive_mass_dt(vars.reproduction_mass_fraction,
                                         vars.net_mass_production);
  vars.dtotal_mass_dt      = s->dtotal_mass_dt(vars.reproduction_mass_fraction,
                                          vars.net_mass_production,
                                          vars.heartwood_mass_rate);
  vars.dabove_ground_mass_dt =
    s->dabove_ground_mass_dt(leaf_area,
                             vars.reproduction_mass_fraction,
                             vars.net_mass_production,
                             vars.heartwood_mass_rate,
                             leaf_area_growth_rate);

  // Odd one out:
  vars.dstem_diameter_dstem_area =
    s->dstem_diameter_dstem_area(vars.bark_area,
                               vars.sapwood_area,
                               vars.heartwood_area);
}

// [eqn 20] Survival of seedlings during germination
double Plant::germination_probability(const Environment& environment) {
  return strategy->germination_probability(environment);
}

// ODE interface -- note that the don't care about time in the plant;
// only Patch and above does.
ode::const_iterator Plant::set_ode_state(ode::const_iterator it) {
  set_height(*it++);
  set_mortality(*it++);
  set_fecundity(*it++);
  set_heartwood_area(*it++);
  set_heartwood_mass(*it++);
  return it;
}
ode::iterator Plant::ode_state(ode::iterator it) const {
  *it++ = height();
  *it++ = mortality();
  *it++ = fecundity();
  *it++ = heartwood_area();
  *it++ = heartwood_mass();
  return it;
}
ode::iterator Plant::ode_rates(ode::iterator it) const {
  *it++ = height_rate();
  *it++ = mortality_rate();
  *it++ = fecundity_rate();
  *it++ = heartwood_area_rate();
  *it++ = heartwood_mass_rate();
  return it;
}

std::vector<std::string> Plant::ode_names() {
  return std::vector<std::string>({"height", "mortality", "fecundity",
	"heartwood_area", "heartwood_mass"});
}

// * R interface
Plant::strategy_type Plant::r_get_strategy() const {
  return *strategy.get();
}

Plant Plant::r_copy() const {
  return *this;
}

Plant_internals Plant::r_internals() const {
  return vars;
}

// This exists only so that I know that nothing will change the
// control parameters by only being able to access a const reference
// (it's shared with everything else that shares the strategy).  It
// also saves a little ugly looking referencing.
const Control& Plant::control() const {
  return strategy->control;
}

// * Private methods

// * Individual size
// [eqn 1-8] Update size variables to a new leaf mass.
void Plant::compute_vars_size(double height_) {
  vars.height = height_;
  vars.leaf_area = strategy->leaf_area(vars.height);
  vars.leaf_mass = strategy->leaf_mass(vars.leaf_area);
  vars.sapwood_area =  strategy->sapwood_area(vars.leaf_area);
  vars.sapwood_mass =  strategy->sapwood_mass(vars.sapwood_area, vars.height);
  vars.bark_area =  strategy->bark_area(vars.leaf_area);
  vars.bark_mass = strategy->bark_mass(vars.bark_area, vars.height);
  vars.root_mass = strategy->root_mass(vars.leaf_area);
  vars.live_mass = strategy->live_mass(vars.leaf_mass, vars.sapwood_mass,
                                       vars.bark_mass , vars.root_mass);
  vars.stem_area = strategy->stem_area(vars.bark_area, vars.sapwood_area,
                                         vars.heartwood_area);
  vars.total_mass = strategy->total_mass(vars.leaf_mass, vars.bark_mass, vars.sapwood_mass,
                                         vars.heartwood_mass, vars.root_mass);
  vars.above_ground_mass = strategy->above_ground_mass(vars.leaf_mass, vars.bark_mass,
                                                       vars.sapwood_mass, vars.root_mass);
  vars.stem_diameter = strategy->stem_diameter(vars.stem_area);
}

Plant make_plant(Plant::strategy_type s) {
  return Plant(make_strategy_ptr(s));
}

}
