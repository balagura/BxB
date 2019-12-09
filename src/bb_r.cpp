//
// B*B simulation of beam-beam effects in van der Meer scans at LHC.
// Copyright (C) 2019 Vladislav Balagura (balagura@cern.ch)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------
//
//  Interface to R language
//
#include <Rcpp.h>
#include "bb.hh"

using namespace Rcpp;
using namespace std;

//' @name BxB
//'
//' @title B*B, beam-beam simulation for vdM scans at LHC
//'
//' @description "beam_beam" function simulates bunch particles in an
//'              accelerator influenced by the electromagnetic force of
//'              another bunch ("beam-beam" effect). Reports corrections
//'              to the bunches overlap integral. See "?beam_beam" for details.
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
//' @param kicked List with the parameters of the kicked bunch. It can conveniently be 
//'           created with the \code{\link{kicked}} function. See \code{?kicked} for details.
//'
//' @param kickers List for the kicker bunches, can be created using the
//'           \code{\link{kickers}} function. See \code{?kickers} for details.
//'
//' @param sim List with B*B simulation options. Can be conveniently created with the
//'           \code{\link{sim}} function. See \code{?sim} for details.
//'
//' @param quiet controls whether the input and the execution progress will be printed to cout
//'
//' @return Results of the simulation as data.frame with the
//' fields\preformatted{"ip"          - interaction point, counting from zero
//' "step"        - scan step number, counting from zero
//' "correction" - beam-beam/no beam-beam luminosity correction
//' "no_bb_analytic_integ"                  - analytic overlap integral
//' "no_bb_numeric_over_analytic_integ"     - numeric/analytic ratio
//' "no_bb_numeric_over_analytic_integ_err" - its error, all without beam-beam
//' "avr_x_analytic", "avr_y_analytic"           - analytic center-of-mass shift
//' "avr_x_numeric", "avr_y_numeric"             - numeric shift
//' "no_bb_avr_x_numeric", "no_bb_avr_y_numeric" - analytic shift without beam-beam}
//'
//' In addition, depending on the options specified in \code{sim.output},
//' the program can print out several files in the \code{sim.output_dir}
//' directory. See \code{\link{?sim}} for details about that and the simulation
//' in general.
//'
//' @author Vladislav BALAGURA <balagura@cern.ch>
//'
//' @examples
//' beam_beam(kicked = kicked(momentum = 3500, Z=1, ip=1,
//'                           beta = list(x=c(rep(1.5,3),6),
//'                                       y=c(rep(1.5,3),6)), 
//'                           next_phase_over_2pi = list(x = 0.31*(1:4),
//'                                                      y = 0.32*(1:4)),
//'                           gaussian = list(x = list(sig=rep(40,2), w=c(0.2, 0.8)),
//'                                           y = list(sig=c(39.99, 40), w=c(0.3, 0.7))),
//'                           exact_phases = FALSE),
//'           kickers =  kickers(Z = 1,
//'                              n_particles = rep(8.5e10, 4),
//'                              gaussian = list(x = list(list(sig=rep(40, 2), w=c(0.2, 0.8)),
//'                                                       list(sig=rep(40, 3), w=c(0.3, 0.6, 0.1)),
//'                                                       list(sig=c(40.002, 40.001, 40.001, 39.998),
//'                                                            w = c(2, 10, 10, 2)),
//'                                                       list(sig=80, w=1)),
//'                                              y = list(list(sig=40, w=0.2),
//'                                                       list(sig=40, w=0.2),
//'                                                       list(sig=40, w=0.2),
//'                                                       list(sig=c(80.001, 79.999), w=rep(1, 2)))),
//'                              position = list(x = list(10*(0:20), 10*(0:20), 10*(0:20), 10*(0:20)),
//'                                              y = list(0,0,0,0))),
//'           sim =  sim(n_points = 5000,
//'                      ns = c(1000, 1000, 0, 5000),
//'                      kick_model = 'precise',
//'                      n_sigma_cut = 5,
//'                      density_and_field_interpolators_n_cells_along_grid_side = c(500, 500),
//'                      n_random_points_to_check_interpolation = 10000,
//'                      select_one_turn_out_of = 1000,
//'                      seed = 123456789,
//'                      output_dir = "",
//'                      output = "")
//'
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
  Sim s;
  s.n_points = as<int>(sim["n_points"]);
  const auto& n_turns = as<std::vector<int> >(sim["n_turns"]);
  if (n_turns.size() != PHASES) {
    Rcpp::stop("The length of \"n_turns\" array must be " + to_string(PHASES));
  }
  copy(n_turns.begin(), n_turns.end(), s.n_turns);
  s.kick_model = as<string>(sim["kick_model"]);
  s.n_sigma_cut = as<double>(sim["n_sigma_cut"]);
  vector<int> n_cells = as<vector<int> >(sim["density_and_field_interpolators_n_cells_along_grid_side"]);
  if (n_cells.size() != 2) {
    Rcpp::stop("The length of \"density_and_field_interpolators_n_cells_along_grid_side\" array must be 2");
  }
  copy(n_cells.begin(), n_cells.end(),
       s.density_and_field_interpolators_n_cells_along_grid_side);
  s.n_random_points_to_check_interpolation = as<int>(sim["n_random_points_to_check_interpolation"]);
  s.select_one_turn_out_of = as<int>(sim["select_one_turn_out_of"]);
  s.seed = as<long int>(sim["seed"]);
  s.output_dir = as<string>(sim["output_dir"]);
  s.output = as<string>(sim["output"]);
  vector<vector<BB_Summary_Per_Step_IP> > sum;
  beam_beam(k, ks, s, &sum, quiet);
  int n_step = sum[0].size();
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
      no_bb_avr_x_numeric[i] = bb.no_bb_avr_numeric[0];
      no_bb_avr_y_numeric[i] = bb.no_bb_avr_numeric[1];
      avr_x_analytic[i] = bb.avr_analytic[0];
      avr_y_analytic[i] = bb.avr_analytic[1];
      avr_x_numeric[i] = bb.avr_numeric[0];
      avr_y_numeric[i] = bb.avr_numeric[1];
    }
  }
  return DataFrame::create(_["ip"] = v_ip,
			   _["step"] = v_step,
			   _["correction"] = correction,
			   _["no_bb_analytic_integ"] = no_bb_analytic_integ,
			   _["no_bb_numeric_over_analytic_integ"] = no_bb_numeric_over_analytic_integ,
			   _["no_bb_numeric_over_analytic_integ_err"] = no_bb_numeric_over_analytic_integ_err,
			   _["no_bb_avr_x_numeric"] = no_bb_avr_x_numeric,
			   _["no_bb_avr_y_numeric"] = no_bb_avr_y_numeric,
			   _["avr_x_analytic"] = avr_x_analytic,
			   _["avr_y_analytic"] = avr_y_analytic,
			   _["avr_x_numeric"] = avr_x_numeric,
			   _["avr_y_numeric"] = avr_y_numeric);
}
