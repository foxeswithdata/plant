// Built from  src/ff16_cohort.cpp on Wed Sep 16 07:11:11 2020 using the scaffolder, from the strategy:  FF16
#include <plant.h>

// Helpers for ES20 model

// [[Rcpp::export]]
double cohort_schedule_max_time_default__Parameters___ES20__ES20_Env(const plant::Parameters<plant::ES20_Strategy,plant::ES20_Environment>& p) {
   return plant::cohort_schedule_max_time_default<plant::Parameters<plant::ES20_Strategy,plant::ES20_Environment> >(p);
}

// [[Rcpp::export]]
plant::CohortSchedule cohort_schedule_default__Parameters___ES20__ES20_Env(const plant::Parameters<plant::ES20_Strategy,plant::ES20_Environment>& p) {
   return plant::cohort_schedule_default<plant::Parameters<plant::ES20_Strategy,plant::ES20_Environment> >(p);
}

// [[Rcpp::export]]
plant::CohortSchedule make_cohort_schedule__Parameters___ES20__ES20_Env(const plant::Parameters<plant::ES20_Strategy,plant::ES20_Environment>& p) {
   return plant::make_cohort_schedule<plant::Parameters<plant::ES20_Strategy,plant::ES20_Environment> >(p);
}
