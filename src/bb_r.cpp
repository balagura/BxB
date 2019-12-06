#include <Rcpp.h>
#include "bb.hh"

using namespace Rcpp;
using namespace std;

//' @name BxB
//'
//' @title B*B, beam-beam simulation for vdM scans at LHC
//'
//' @description "beam_beam" function simulates bunch particles in an
//'              accelerator influenced by the electromagnetic interaction
//'              with another bunch ("beam-beam" effect). Reports corrections
//'              to the bunches overlap integral.
//'
//' @docType package
//' @name BxB
//' @author Vladislav BALAGURA <balagura@cern.ch>
//'
//' @import Rcpp
//' @useDynLib BxB
	    
void gaussian(SEXP rg, vector<Gaussian>& g) {
  vector<double>
    sig = as<std::vector<double> >((as<List>(rg))["sig"]),
    w   = as<std::vector<double> >((as<List>(rg))["w"]);
  if (sig.size() != w.size()) {
    Rcpp::stop("given numbers of Gaussian sigmas and weights differ");
  }
  g.resize(sig.size());
  for (size_t gauss=0; gauss<sig.size(); ++gauss) {
    g[gauss].sig = sig[gauss];
    g[gauss].w   = w  [gauss];
  }
}

//' @title beam_beam
//'
//' @description Simulates bunch particles in an accelerator influenced by the
//'              electromagnetic interaction with another bunch ("beam-beam"
//'              effect). Reports corrections to the bunches overlap integral.
//'
//' @return data.frame summary, including corrections to overlap integrals
//'
//' @author Vladislav BALAGURA <balagura@cern.ch>
//' @export
// [[Rcpp::export]]
List beam_beam(List kicked, List kickers, List sim, bool quiet = false) {
  Kicked k;
  k.momentum = as<double>(kicked["momentum"]);
  k.Z = as<int>(kicked["Z"]);
  k.ip = as<int>(kicked["ip"]);
  string xy_str[2] = {"x", "y"};
  for (int coor=0; coor<2; ++coor) {
    const char* xy = xy_str[coor].c_str();
    k.beta               [coor] = as<vector<double> >(as<List>(kicked["beta"])[xy]);
    k.next_phase_over_2pi[coor] = as<vector<double> >(as<List>(kicked["next_phase_over_2pi"])[xy]);
    gaussian(as<List>(kicked["gaussian"])[xy], k.gaussian[coor]);
  }
  k.exact_phases = as<bool>(kicked["exact_phases"]);
  //
  size_t n_ip = k.beta[0].size();
  Kickers ks;
  ks.Z = as<int>(kickers["Z"]);
  ks.n_particles = as<vector<double> >(kickers["n_particles"]);
  for (int coor=0; coor<2; ++coor) {
    const char* xy = xy_str[coor].c_str();
    ks.gaussian[coor].resize(n_ip);
    ks.position[coor].resize(n_ip);
    for (size_t ip=0; ip<n_ip; ++ip) {
      gaussian(as<List>(as<List>(kickers["gaussian"])[xy])[ip], ks.gaussian[coor][ip]);
      ks.position[coor][ip] = as<vector<double> >(as<List>(as<List>(kickers["position"])[xy])[ip]);
    }
  }
  size_t n_step = ks.position[0][0].size();
  Sim s;
  s.n_points = as<int>(sim["n_points"]);
  const auto& nturns = as<std::vector<int> >(sim["n_turns"]);
  copy(nturns.begin(), nturns.end(), s.n_turns);
  s.kick_model = as<string>(sim["kick_model"]);
  s.n_sigma_cut = as<double>(sim["n_sigma_cut"]);
  vector<int> n_cells = as<vector<int> >(sim["density_and_field_interpolators_n_cells_along_grid_side"]);
  s.density_and_field_interpolators_n_cells_along_grid_side[0] = n_cells[0];
  s.density_and_field_interpolators_n_cells_along_grid_side[1] = n_cells[1];
  s.n_random_points_to_check_interpolation = as<int>(sim["n_random_points_to_check_interpolation"]);
  s.select_one_turn_out_of = as<int>(sim["select_one_turn_out_of"]);
  s.seed = as<long int>(sim["seed"]);
  s.output_dir = as<string>(sim["output_dir"]);
  s.output = as<string>(sim["output"]);
  vector<vector<BB_Summary_Per_Step_IP> > sum(n_ip, vector<BB_Summary_Per_Step_IP>(n_step));
  beam_beam(k, ks, s, &sum, quiet);
  vector<double>
    v_ip(n_step * n_ip),
    v_step(n_step * n_ip),
    correction(n_step * n_ip),
    no_bb_analytic_integ(n_step * n_ip),
    no_bb_numeric_over_analytic_integ(n_step * n_ip),
    no_bb_numeric_over_analytic_integ_err(n_step * n_ip),
    avr_x_analytic(n_step * n_ip),
    avr_y_analytic(n_step * n_ip),
    avr_x_numeric(n_step * n_ip),
    avr_y_numeric(n_step * n_ip),
    no_bb_avr_x_numeric(n_step * n_ip),
    no_bb_avr_y_numeric(n_step * n_ip);
  for (size_t ip=0; ip<n_ip; ++ip) {
    for (size_t step=0; step<n_step; ++step) {
      const BB_Summary_Per_Step_IP& bb = sum[ip][step];
      int i = ip * n_step + step;
      v_ip[i] = ip;
      v_step[i] = step;
      correction[i] = bb.correction;
      no_bb_analytic_integ[i] = bb.no_bb_analytic_integ;
      no_bb_numeric_over_analytic_integ[i] = bb.no_bb_numeric_over_analytic_integ;
      no_bb_numeric_over_analytic_integ_err[i] = bb.no_bb_numeric_over_analytic_integ_err;
      avr_x_analytic[i] = bb.avr_analytic[0];
      avr_y_analytic[i] = bb.avr_analytic[1];
      avr_x_numeric[i] = bb.avr_numeric[0];
      avr_y_numeric[i] = bb.avr_numeric[1];
      no_bb_avr_x_numeric[i] = bb.no_bb_avr_numeric[0];
      no_bb_avr_y_numeric[i] = bb.no_bb_avr_numeric[1];
    }
  }
  return DataFrame::create(_["ip"] = v_ip,
			   _["step"] = v_step,
			   _["correction"] = correction,
			   _["no_bb_analytic_integ"] = no_bb_analytic_integ,
			   _["no_bb_numeric_over_analytic_integ"] = no_bb_numeric_over_analytic_integ,
			   _["no_bb_numeric_over_analytic_integ_err"] = no_bb_numeric_over_analytic_integ_err,
			   _["avr_x_analytic"] = avr_x_analytic,
			   _["avr_y_analytic"] = avr_y_analytic,
			   _["avr_x_numeric"] = avr_x_numeric,
			   _["avr_y_numeric"] = avr_y_numeric,
			   _["no_bb_avr_x_numeric"] = no_bb_avr_x_numeric,
			   _["no_bb_avr_y_numeric"] = no_bb_avr_y_numeric);
}
