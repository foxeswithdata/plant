#include "metacommunity.h"

namespace model {

SEXP metacommunity(Rcpp::CppClass individual, Parameters p) {
  std::string individual_type =
    util::rcpp_class_demangle(Rcpp::as<std::string>(individual));
  SEXP ret = R_NilValue;
  if (individual_type == "Plant") {
    Metacommunity<Plant> obj(p);
    ret = Rcpp::wrap(obj);
  } else if (individual_type == "CohortDiscrete") {
    Metacommunity<CohortDiscrete> obj(p);
    ret = Rcpp::wrap(obj);
  } else {
    ::Rf_error("Cannot make Metacommunity of %s",
	       individual_type.c_str());
  }
  return ret;
}

}