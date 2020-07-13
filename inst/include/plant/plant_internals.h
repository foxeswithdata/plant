// -*-c++-*-
#ifndef PLANT_PLANT_PLANT_INTERNALS_H_
#define PLANT_PLANT_PLANT_INTERNALS_H_

namespace plant {

// These are common to all minimal plants, for now at least.
//
// Moving to a more general "size" based model would be easy enough
// but we'd need to also store height because Patch & Environment
// between them use height to work out how far up to compute the
// canopy openness for.  So like leaf_area being carried around we'd
// need to carry height as well.
struct Plant_internals {
  Plant_internals()
    :
    height(NA_REAL),
    height_dt(NA_REAL),
    mortality(0.0),
    mortality_dt(NA_REAL),
    fecundity(0.0),
    fecundity_dt(NA_REAL),
    area_heartwood(0.0),
    area_heartwood_dt(NA_REAL),
    mass_heartwood(0.0),
    mass_heartwood_dt(NA_REAL),
    area_leaf(NA_REAL),
    area_leaf_dt(0.0)
    mass_leaf(NA_REAL),
    mass_leaf_dt(0.0),
    area_sapwood(NA_REAL),
    area_sapwood_dt(0.0)
    mass_sapwood(NA_REAL),
    mass_sapwood_dt(0.0),
    area_bark(NA_REAL),
    area_bark_dt(0.0),
    mass_bark(NA_REAL),
    mass_bark_dt(0.0),
    mass_root(NA_REAL),
    mass_root_dt(0.0),
    mass_store(NA_REAL),
    mass_store_dt(0.0),
    area_stem(NA_REAL),
    area_stem_dt(0.0)
    {
  }
  double height;
  double area_leaf;
  double area_leaf_dt;
  double height_dt;
  double mortality;
  double mortality_dt;
  double fecundity;
  double fecundity_dt;
  double area_heartwood;
  double area_heartwood_dt;
  double mass_heartwood;
  double mass_heartwood_dt;
  double area_sapwood;
  double area_sapwood_dt;
  double mass_sapwood;
  double mass_sapwood_dt;
  double area_bark;
  double area_bark_dt;
  double mass_bark;
  double mass_bark_dt;  
  double mass_root;
  double mass_root_dt;  
  double mass_store;
  double mass_store_dt;  
  double area_stem;
  double area_stem_dt;
  
};

}

#endif
