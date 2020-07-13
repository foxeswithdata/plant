#include <plant/environment.h>
#include <plant/parameters.h>

#include <random>

namespace plant {

Environment::Environment(double disturbance_mean_interval,
                         std::vector<double> seed_rain_,
                         Control control)
  : time(0.0),
    disturbance_regime(disturbance_mean_interval),
    seed_rain(seed_rain_),
    seed_rain_index(0),
    light_environment_generator(make_interpolator(control)) {

    stress_mean = 0.684931506849315;
    stress_sd = 0.0684931506849315;
    prepare_stress();
    test = 0.8;
}

double Environment::canopy_openness(double height) const {
  const bool within_canopy = height <= light_environment.max();
  return within_canopy ? light_environment.eval(height) : 1.0;
}

// Computes the probability of survival from 0 to time.
double Environment::patch_survival() const {
  return disturbance_regime.pr_survival(time);
}

// Computes the probability of survival from time_at_birth to time, by
// conditioning survival over [0,time] on survival over
// [0,time_at_birth].
double Environment::patch_survival_conditional(double time_at_birth) const {
  return disturbance_regime.pr_survival_conditional(time, time_at_birth);
}

// Reset the environment.
void Environment::clear() {
  time = 0.0;
  clear_light_environment();
}

void Environment::clear_light_environment() {
  light_environment.clear();
}

double Environment::seed_rain_dt() const {
  if (seed_rain.empty()) {
    Rcpp::stop("Cannot get seed rain for empty environment");
  }
  return seed_rain[seed_rain_index];
}

void Environment::set_seed_rain_index(size_t x) {
  seed_rain_index = x;
}




// stress interface 
double Environment::time_in_year() const{
  return (time-floor(time));
}

bool Environment::stressed() const {
  double yr = floor(time);
  if(time_in_year() < stress_regime[yr]){
    return false;
  }
  else{
    return true;
  }
}

void Environment::prepare_stress(){
  std::default_random_engine stress_regime_engine;
  
  // values near the mean are the most likely
  // standard deviation affects the dispersion of generated values from the mean
  std::normal_distribution<double> d(stress_mean,stress_sd);
  
  //hard coded 3000 as an upper limit that will probably not be reached
  for (int i = 0; i < 3000; i++) {
    stress_regime[i] = d(stress_regime_engine);
  }
}

// * R interface
void Environment::r_set_seed_rain_index(util::index x) {
  set_seed_rain_index(x.check_bounds(seed_rain.size()));
}

void prepare_stress(){
  std::default_random_engine stress_regime_engine;
  
  // values near the mean are the most likely
  // standard deviation affects the dispersion of generated values from the mean
  std::normal_distribution<double> d(stress_mean,stress_sd);
  
  //hard coded 3000 as an upper limit that will probably not be reached
  for (int i = 0; i < 3000; i++) {
    stress_regime[i] = d(stress_regime_engine);
  }
}

}
