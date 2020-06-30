// #include <plant/uniroot.h>
// #include <plant/qag.h>
// #include <plant/models/assimilation.h>
// #include <plant/models/NEW_strategy.h> // Reference your file
// #include <RcppCommon.h> // NA_REAL
// 
// namespace plant {
// 
//     // Constructor
//     NEW_Strategy::NEW_Strategy() {
//       // initialise parameters
//       param_1 = -9999;
//       param_2 = -9999;
//       param_3 = -9999;
// 
//       collect_all_auxillary = false;
// 
//       refresh_indices();
//       name = "strategy_name";
//     }
// 
//       /***
//        * No need to touch this
//        */
//       void NEW_Strategy::refresh_indices () {
//       // Create and fill the name to state index maps
//       state_index = std::map<std::string,int>();
//       aux_index   = std::map<std::string,int>();
//       std::vector<std::string> aux_names_vec = aux_names();
//       std::vector<std::string> state_names_vec = state_names();
//       for (int i = 0; i < state_names_vec.size(); i++) {
//       state_index[state_names_vec[i]] = i;
//       }
//       for (int i = 0; i < aux_names_vec.size(); i++) {
//       aux_index[aux_names_vec[i]] = i;
//       }
//       }
// 
//       /***
//        * Implement the model
//        */
//       double NEW_STRATEGY::equation_1(double arg_1, arg_2) {
//         double var = -9999;
//         return var;
//       }
// 
//       /***
//        * Updating dependent auxilary states
//        * Need to figure this out further
//        */
//       void NEW_Strategy::update_dependent_aux(const int index, Internals& vars) {
//       if (index == HEIGHT_INDEX) {
//       double height = vars.state(HEIGHT_INDEX);
//       vars.set_aux(aux_index.at("competition_effect"), area_leaf(height));
//       }
//       }
// 
// 
//       /***
//        * Here is where you call and update the strategy variables
//        *
//        * one-shot update of the scm variables
//        * i.e. setting rates of ode vars from the state and updating aux vars
//        */
//       void NEW_Strategy::compute_rates(const NEW_Environment& environment,
//       bool reuse_intervals,
//       Internals& vars) {
// 
//       // Recover necessary state variables
//       double state_1 = vars.state(STATE_1_INDEX);
//       double state_2 = vars.state(STATE_2_INDEX);
//       double state_3 = vars.state(STATE_3_INDEX);
// 
// 
//       // Recover necessary auxilary state variables
//       double aux_state_1 = vars.aux(aux_index.at("aux_state_1"));
// 
//       const double aux_state_2_dt_ = equation_1();
// 
//       // store the aux sate
//       vars.set_aux(aux_index.at("aux_state_2_dt"), aux_state_2_dt_);
// 
//       //
//       // Work through the model
//       //
// 
//       // set new rates
//       vars.set_rate(STATE_1_INDEX, 0.0);
//       vars.set_rate(STATE_2_INDEX, 0.0);
//       vars.set_rate(STATE_3_INDEX, 0.0);
// 
//       }
// 
//     /***
//      * Set up the strategy
//      */
//     void NEW_Strategy::prepare_strategy() {
//     // Set up the integrator
//     control.initialize();
//     assimilator.initialize(a_p1, a_p2, eta);
// 
//     // pre-compute any variable that's needed
// 
//     }
// 
//     /***
//      * Leave as is
//      */
//     NEW_Strategy::ptr make_strategy_ptr(NEW_Strategy s) {
//     s.prepare_strategy();
//     return std::make_shared<NEW_Strategy>(s);
//     }
// }
//   
