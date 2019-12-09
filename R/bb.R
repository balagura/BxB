#
# B*B simulation of beam-beam effects in van der Meer scans at LHC.
# Copyright (C) 2019 Vladislav Balagura (balagura@cern.ch)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------
#
# Convenience functions to form kicked, kickers and sim arguments for
# beam_beam(kicked, kickers, sim) B*B simulation function.
#

#' @title \code{kicked} convenience function
#' 
#' @description Convenience function aggregating parameters of the kicked
#'     bunch into the "kicked" argument of the
#'     \code{\link{beam_beam}(kicked, kickers, sim, quiet)}" simulation function.
#'
#' @param momentum In GeV.
#' 
#' @param Z Charge of the particles in the units of the proton charge.
#' 
#' @param ip Interaction point (IP) number (counting from zero) for which the
#'     Gaussian sigmas are specified. Kicked bunch sigmas at other IP2 are
#'     calculated via beta-function values as sigma(ip) * sqrt(beta(ip) /
#'     beta(IP2). Note, that contrary to that, the kicker sigmas should be
#'     specified directly at the IP of the beam-beam interaction, without any
#'     extrapolation via beta. In other words, such kicker sigmas, in general,
#'     can differ from the ones at this kicked "ip".
#'
#' @param beta List with "x" and "y" components each specifying a vector of
#'     the corresponding beta-function values (beta-star) at all simulated
#'     interaction points, in meters.
#' 
#' @param next_phase_over_2pi List with "x" and "y" components specifying
#'     transverse betatron oscillation phases divided by 2pi.  The phase at
#'     IP=0 is zero by definition, therefore, the "x" and "y" vectors should
#'     start from the phase at IP=1, this is why this parameter is called
#'     "next" phase over 2pi. The length of "x" and "y" vectors should be
#'     equal to the total number of simulated IPs.
#'
#' @param gaussian List with "x" and "y" components. Each component is in turn
#'     a list with two elements, "sig" and "w" specifying Gaussian sig(mas) in
#'     um and the corresponding w(eights) of the multi-Gaussian kicked bunch
#'     density. All Gaussians should have a common mean at zero. The weights
#'     might be given not normalized.
#' 
#' @param exact_phases If FALSE (default), small irrational numbers will be
#'     added to all phases to increase the randomness in the transverse
#'     particle trajectories and improve the simulation, see \code{Details}.
#'
#' @details The returned list contains all information on the kicked bunch
#'     needed for the B*B simulation.
#' 
#' If values in "next_phase_over_2pi" are given with 2-3 digits after the
#' comma, after 100 or 1000 turns the points return almost to their original
#' positions if the beam-beam effect is small (as 100* or 1000 * phases/2pi
#' become integer). The simulation, therefore, probes the same part of the
#' phase space over and over again, and extra turns almost do not improve the
#' convergence. To avoid this and to add extra randomness to the simulation,
#' "negligible" irrational numbers randomly distributed in the interval
#' -exp(-8.5)...+exp(-8.5) = +/-0.0001017342 are added by default to all phase
#' / 2pi values. This makes them irrational and opens otherwise closed
#' Lissajous figures of betatron oscillations in X-Y plane. In reality the
#' phases/2pi are always irrational.
#'
#' If this is not desired, an option "exact.phases = TRUE" can be set. Then
#' all phases/2pi are used exactly as they are given.
#' 
#' @return List with all function arguments as list components. It contains
#'     full information on the kicked bunch required for the B*B
#'     simulation. The list can be fed to \code{\link{beam_beam}(kicked,
#'     kickers, sim, quiet)} simulation function.
#'
#' @seealso \code{kickers}, \code{sim}, \code{beam_beam}.
#'
#'@examples
#' kicked(momentum = 3500, Z=1, ip=1,
#'        beta = list(x=c(rep(1.5,3),6),
#'                    y=c(rep(1.5,3),6)), 
#'        next_phase_over_2pi = list(x = 0.31*(1:4),
#'                                   y = 0.32*(1:4)),
#'        gaussian = list(x = list(sig=rep(40,2), w=c(0.2, 0.8)),
#'                        y = list(sig=c(39.99, 40), w=c(0.3, 0.7))),
#'        exact_phases = FALSE)
#'
#' 
#' @author Vladislav BALAGURA <balagura@cern.ch>
#' @export
# [[Rcpp::export]]
kicked <- function(momentum, Z, ip, beta, next_phase_over_2pi, gaussian, exact_phases = FALSE) {
    list(momentum = momentum, Z = Z, ip = ip,
         beta = beta, next_phase_over_2pi = next_phase_over_2pi,
         gaussian = gaussian, exact_phases = exact_phases)
}

#' @title \code{kickers} convenience function
#' 
#' @description Convenience function aggregating parameters of the kicker
#'     bunches into the "kickers" argument of the
#'     \code{\link{beam_beam}(kicked, kickers, sim, quiet)}" simulation function.
#'
#' @param Z Charge of the particles in the units of the proton charge.
#' 
#' @param n_particles Numeric vector specifying the number of particles in all
#'     kickers
#'
#' @param gaussian List of lists of lists. The outermost has "x" and "y"
#'     components. The number of components in the middle list is equal to the
#'     number of interaction points (IPs). The innermost list has two
#'     elements, "sig" and "w" specifying Gaussian sig(mas) in um and the
#'     corresponding w(eights) of the multi-Gaussian kicker bunch densities at
#'     the given IP. All Gaussians should have a common mean. The weights
#'     might be given not normalized. For example, \preformatted{gaussian = list(x = list(list(sig = 50, w=1),
#'                          list(sig=60, w=10)),
#'                 y = list(list(sig = c(40,60), w=c(1,2)),
#'                          list(sig=50, w=1)))}
#'
#'     specifies kicker densities for the simulation with two IPs.  Note, that
#'     contrary to the kicked bunch, the kicker sigmas should be specified
#'     directly at the IP of the beam-beam interaction, without any
#'     extrapolation via beta. In other words, such kicker sigmas, in general,
#'     can differ from the ones at the \code{ip} where the kicked sigmas are
#'     specified
#'
#' 
#' @param position List of lists. The outer has "x" and "y" components
#'     specifying kicker positions in um in the frame where the bunch center
#'     of the kicked bunch before beam-beam was at zero. The inner list
#'     contains "n_ip" vectors where "n_ip" is the number of IPs. Each vector
#'     can have either one or "n_step" kicker positions for a given IP, where
#'     "n_step" is the number of scan points in van der Meer scan. If the
#'     vector has only one position, it is "recycled", and the kicker position
#'     is assumed to be constant during the scan. For example,\preformatted{position = list(x = list(10*(0:20), 500, 10*(0:20), 500),
#'                 y = list(0,0,0,0))}
#' 
#'     specifies the simultaneous scan at IP=0 and IP=2 along X with the beams
#'     at IP=1, 3 separated by 500 um, and without Y-separations in all IPs.
#' 
#' @details 
#' 
#' @return List with all function arguments as list components. It contains
#'     full information on the kicker bunches required for the B*B
#'     simulation. The list can be fed to \code{\link{beam_beam}(kicked,
#'     kickers, sim, quiet)} simulation function.
#'
#' @seealso \code{kicked}, \code{sim}, \code{beam_beam}.
#'
#' @examples
#' kickers(Z = 1,
#'         n_particles = rep(8.5e10, 4),
#'         gaussian = list(x = list(list(sig=rep(40, 2), w=c(0.2, 0.8)),
#'                                  list(sig=rep(40, 3), w=c(0.3, 0.6, 0.1)),
#'                                  list(sig=c(40.002, 40.001, 40.001, 39.998),
#'                                       w = c(2, 10, 10, 2)),
#'                                  list(sig=80, w=1)),
#'                         y = list(list(sig=40, w=0.2),
#'                                  list(sig=40, w=0.2),
#'                                  list(sig=40, w=0.2),
#'                                  list(sig=c(80.001, 79.999), w=rep(1, 2)))),
#'         position = list(x = list(10*(0:20), 10*(0:20), 10*(0:20), 10*(0:20)),
#'                         y = list(0,0,0,0)))
#' 
#' @author Vladislav BALAGURA <balagura@cern.ch>
#' @export
# [[Rcpp::export]]
kickers <- function(Z, n_particles, gaussian, position) {
    list(Z = Z, n_particles = n_particles, gaussian = gaussian, position = position)
}

#' @title \code{sim} convenience function
#' 
#' @description Convenience function aggregating parameters controlling B*B
#'     simulation into the "sim" argument of the
#'     \code{\link{beam_beam}(kicked, kickers, sim, quiet)}" function.
#'
#' @param n_points Number of "macro-particles" traced in the simulation will
#'     be slightly less than this number (eg. by 21\% for round single Gaussian
#'     kicked bunch). See \code{Initial distribution of macro-particles}.
#' 
#' @param n_turns Vector with 4 integers: number of turns in the
#'     "no beam-beam", "adiabatic switch on", "stabilization" and "beam-beam"
#'     simulation phases.
#' 
#' @param kick_model One of "precise", "average" or
#'     "precise.minus.average". If "precise" (default), the kick is calculated
#'     according to the exact formula. "average" sets the kick to the constant
#'     not depending on X,Y-position of the particle. This is equivalent to
#'     the field of a dipole magnet. Its strength is chosen to reproduce the
#'     precise overall force on the kicked bunch. "precise.minus.average" sets
#'     the kick to the difference between the "precise" and "average" values.
#'
#' @param n_sigma_cut Controls the limits of 4D phase space volume from which
#'     the macro-particles are selected. See \code{Initial distribution of
#'     macro-particles}.
#' 
#' @param density_and_field_interpolators_n_cells_along_grid_side Two integers
#'     specifying the number of cells in the interpolation grids. The first is
#'     for a one-dimensional interpolation of the kicker bunch densities in X
#'     and Y, the second is for each side of the two-dimensional X-Y
#'     interpolation grid of the kicker fields. If one of these numbers is set
#'     to zero, the corresponding interpolation is not performed and the value
#'     each time is calculated exactly (CPU consuming). See
#'     \code{Interpolation}.
#' 
#' @param n_random_points_to_check_interpolation If this parameter is larger
#'     than zero and \code{quiet} parameter in \code{\link{beam_beam}}
#'     simulation function is not \code{TRUE}, the quality of the
#'     interpolations will be estimated by comparing the mismatches between
#'     the exact and the interpolated values at
#'     \code{n_random_points_to_check_interpolation} points randomly
#'     distributed inside the interpolation grid. The results will be printed
#'     to \code{cout}.
#' 
#' @param select_one_turn_out_of If \code{output} option \code{points} is
#'     specified, all macro-particle positions will be stored to disk. Storage
#'     per accelerator turn might require too much space. Instead, one can
#'     \code{select_one_turn_out_of} N turns.  Eg. if this parameter is set to
#'     1000, only the positions before the turns = 999, 1999, 2999, ...
#'     (counting from zero) will be stored. If \code{points} option is not
#'     requested in \code{output}, \code{select_one_turn_out_of} has no
#'     effect.
#'
#' @param seed Random seed for the simulation. If specified, the results of
#'     the simulation will be fully reproducible. Otherwise, the current time
#'     will be used as a seed, and the results will be not reproducible.
#'
#' 
#' @param output_dir The name of the subdirectory (relative to the current
#'     path or absolute) where all simulated files will be stored. This
#'     subdirectory is created by the program and should not exist, otherwise
#'     the program will terminate without overwriting anything.
#'
#'
#' @param output A character string with white-space separated options
#'     controlling what will be calculated during the simulation. If
#'     \code{output_dir} is not empty (""), for every option \code{XXXX.XXXX},
#'     one compressed file will be created under the name
#'     \code{output_dir/XXXX_XXXX.txt.gz}. \code{output} might include any
#'     combination of \code{integrals.per.turn}, \code{avr.xy.per.turn},
#'     \code{integrals.per.particle}, \code{avr.xy.per.particle},
#'     \code{points} or be empty.
#' 
#' @section Initial distribution of macro-particles:
#' The traced macro-particles of the "kicked" bunch are selected in the
#' following way.
#'
#' It is assumed that the bunch densities factorize in X and Y. Their
#' projections to X and Y are multi-Gaussian distributions with common
#' centers. Because of the circular motion in X-X' and Y-Y' planes (where
#' X',Y' denote the angular coordinates scaled by the accelerator
#' beta-function, eg. X' = dX/dZ * beta), the projections to X' (Y') are
#' identical to X (Y). So, eg. a single Gaussian in X makes a two-dimensional
#' Gaussian with equal sigmas in X-X' plane, and a multi-Gaussian makes a
#' multi two-dimensional Gaussian with the same weights.
#'
#' The four-dimensional kicked bunch density is sampled internally in a
#' two-dimensional rX-rY grid of the radii in the X-X', Y-Y' planes at IP=0.
#'
#' Each grid side (rX or rY) extends from 0 to a maximal radius rX,Y_max. The
#' latter is chosen such that the circle with the radius rX_max (or rY_max)
#' contains the same fraction of a multi two-dimensional Gaussian as that of a
#' single two-dimensional Gaussian inside the circle with the radius
#' \code{n_sigma_cut}.
#'
#' The number of grid lines in rX and rY are chosen to be the same and equal to
#' int(sqrt(\code{n_points})). Finally, only rX-rY points inside the ellipse
#' inscribed to the rectangle (rX_max X rY_max) are used for the
#' simulation. For every selected rX-rY pair, one point at the X-X', Y-Y'
#' circles with these radii and random phases is traced in the accelerator.
#'
#' @section Four phases in the simulation:
#' First phase is without beam-beam. It allows to calculate numerically the
#' undisturbed overlap integral, compare it with the exact analytic formula and
#' estimate the bias of the numerical integration. The final correction is then
#' calculated as a ratio of the numerical integrals with (last phase) and
#' without (first phase) the beam-beam interaction. The bias is cancelled in
#' the ratio at least partially and this potentially allows to improve the
#' precision.
#'
#' During the second phase the beam-beam interaction is switched on
#' "adiabatically": linearly with the turn number from zero to its nominal
#' value. If this phase is omitted (zero turns) the switch is abrupt: from
#' zero to nominal.
#'
#' After beam-beam is fully switched on, there might be a need to wait for the
#' stabilization. This is the third phase. With the adiabatic switch on
#' (eg. during 1000 turns) the stabilization phase can normally be omitted and
#' the corresponding number of turns can be set to zero.
#'
#' The final, fourth phase is with beam-beam. Only this phase is used to
#' calculate the perturbed overlap integral and the luminosity correction.
#'
#' \code{n_turns} parameter specifies the number of turns in all 4 phases.
#' 
#' @section Interpolation:
#' To speed up the simulation, the X- and Y-densities of the kicker bunch at a
#' given point are determined using linear (one-dimensional) interpolations,
#' while the generated E-field - using a bilinear (two-dimensional)
#' interpolation. For that, the precise density and the field are precalculated
#' at the grid of points and interpolated in between.
#'
#' The interpolation grid is chosen to minimally cover the ranges +/- 1.1 *
#' rX_max, +/- 1.1 * rY_max around the kicked bunch center but viewed from the
#' kicker bunch center for all set kicker separations. Here, rX,Y_max are the
#' maximal simulated radii in X-X' and Y-Y' planes.  Therefore, the grid is
#' sufficiently wide to cover all simulated without beam-beam particle
#' trajectories (circles), while with 10\% margins (with the factor 1.1) it
#' likely covers also the trajectories deformed by the beam-beam
#' interaction. If not, ie. if some point goes beyond the interpolation grid,
#' the program falls back to the calculation of the exact field instead of the
#' interpolation.
#'
#' The linear (bilinear) grid is a raw (a matrix) of N (NxN) cells and has the
#' total number of points N+1 ((N+1)*(N+1)). N is defined separately for the
#' one-dimensional (X,Y)-density and for the two-dimensional field
#' interpolators by the parameter
#' \code{density_and_field_interpolators_n_cells_along_grid_side}. If the
#' interpolation is not desired, the corresponding N should be set to zero,
#' like "0 500". The first number is for the density, the second - for the
#' field, so in this case the density will be calculated using exact formulas
#' while the field will be interpolated with N=500.
#'
#' If the kicker densities and the fields are interpolated, the simulation of
#' complex multi-Gaussian elliptical bunches and simple round
#' single-Gaussian ones takes about the same time.
#'
#' The interpolators are created not for each kicker bunch but for each group
#' of identical bunches. Ie. if there are identical kickers, the interpolators
#' are reused to save memory.
#'
#' The accuracy of the interpolation is estimated if the parameter
#' \code{n_random_points_to_check_interpolation} is specified. This is
#' performed by measuring the maximal and average absolute mismatches between
#' the interpolated and the exact values at
#' \code{n_random_points_to_check_interpolation} distributed inside the
#' interpolation grid. The mismatches are then normalized to the maximal
#' absolute exact value and printed out.
#'
#' If the interpolation was switched off by setting one or both
#' \code{density_and_field_interpolators_n_cells_along_grid_side} values to
#' zero, the corresponding \code{n_random_points_to_check_interpolation} has
#' no effect. The check is not performed either if
#' \code{n_random_points_to_check_interpolation} is set to zero.
#'
#' @section Kick model:
#' For debugging purposes one can redefine the kick formula by setting
#' "kick.model" parameter to one of the following:
#' \preformatted{
#'      precise
#'      precise.minus.average
#'      average
#' }
#' \code{precise} is the default. In this case the exact kick formula is used.
#'
#' \code{average} means constant X,Y-independent kick equal to the value
#' averaged over the kicked bunch. Ie. it is equal to the sum of the kicks of
#' all bunch particles divided by their number. For the Gaussian bunches it
#' can be calculated as the kick exerted by the kicker bunch with (Capital)
#' sigma = sqrt(sigma1^2 + sigma2^2) on one particle placed at the kicked
#' bunch center. The constant kick only shifts the kicked bunch as a whole,
#' but does not modify its shape, so the resulting luminosity change can be
#' computed using analytic formula.
#'
#' \code{precise.minus.average} kick is simply calculated as the difference
#' "precise" - "average". This model allows to compare the "precise" beam-beam
#' luminosity correction with the sum of the corrections obtained with
#' "precise.minus.average" and "average".
#'
#' @section Output: \code{output} parameter controls what should be calculated
#'     and printed after the simulation. It is a white-space separated list of
#'     options listed below:
#'
#'   \code{integrals.per.turn} - overlap integrals per turn (averaged over
#'                        particles). Format: step ip ip.kicker.x ip.kicker.y
#'                        i.turn integral.
#'
#'   \code{avr.xy.per.turn} - average kicked bunch center per turn (averaged
#'                     over particles). Format: step ip i.turn average.x
#'                     average.y.
#'
#'   \code{integrals.per.particle} - overlap integrals per particle, per phase
#'                            (averaged over turns in the given phase),
#'                            assuming that the kicked bunch is squashed to
#'                            this particle, ie. the bunch density is a
#'                            delta-function at its position. Format: step ip
#'                            phase i.particle integral.
#'
#'   \code{avr.xy.per.particle} -   average kicked bunch center per particle, per phase
#'                           as above. Format: step ip phase i.particle avr.x avr.y.
#'
#'   \code{points} -  the traced positions of the kicked bunch particles before
#'             accelerator "i.turn", note, the generated file might be very long.
#'             Format: step ip i.turn i.particle X X' Y Y'.
#'
#' The integers i.turn and i.particle are counted from zero.
#'
#' If \code{output_dir} is not empty (""), for every option \code{XXXX.XXXX},
#' one compressed file will be created under the name
#' \code{output_dir/XXXX_XXXX.txt.gz}. In addition, the files
#' \code{output_dir/rx.ry.weights}, \code{output_dir/kicker_positions.txt} and
#' \code{output_dir/summary.txt} are printed out regardless of options in
#' \code{output} with the following content:
#'
#' \code{output_dir/rx.ry.weights} - for every simulated particle: its initial
#'                   rX and rY radii (in X,X' and Y,Y' phase spaces) at ip and
#'                   its "weight". Format: i.particle ip rx ry weight.
#'
#' \code{output_dir/kicker_positions.txt} - coordinates of the kicker bunch
#'     centers with respect to the nominal (not perturbed by the
#'     beam-beam effect) kicked bunch center. Format: step ip x y.
#'
#'  \code{output_dir/summary.txt} - the main results of the simulation
#'     including the overlap integrals and the corresponding beam-beam
#'     corrections. Format: \preformatted{<step> <IP> <beam-beam/no beam-beam luminosity correction>
#' no beam-beam: <analytic> overlap, <numeric/analytic ratio> and <its error>
#' <X>, <Y> numeric kicked center-of-mass shift without beam-beam
#' <X>, <Y> numeric shift with beam-beam
#' <X>, <Y> analytic shift with beam-beam.}
#'
#' The error of the numeric/analytic overlap ratio without beam-beam is roughly
#' estimated from the turn-by-turn variations, available only if
#' \code{integrals.per.turn} option is set in \code{output}. Otherwise this
#' error is assigned to \code{nan}. Similarly, numeric \code{<X>, <Y>} average
#' shifts are calculated only if \code{avr.xy.per.particle} option is chosen,
#' and set to \code{nan} otherwise. Without beam-beam these shifts should be
#' close to zero.
#'
#' If, minimally, only this summary is required from the simulation, it is
#' better to set \code{output} to am empty string (""). In this case the
#' simulation will run faster, as every option requires extra CPU time.
#' 
#' @return List with all function arguments as list components. The list can
#'     be fed to \code{\link{beam_beam}(kicked, kickers, sim, quiet)}
#'     B*B simulation function.
#' 
#' @seealso \code{kicked}, \code{kickers}, \code{beam_beam}.
#'
#' @examples
#' sim(n_points = 5000,
#'     ns = c(1000, 1000, 0, 5000),
#'     kick_model = 'precise',
#'     n_sigma_cut = 5,
#'     density_and_field_interpolators_n_cells_along_grid_side = c(500, 500),
#'     n_random_points_to_check_interpolation = 10000,
#'     select_one_turn_out_of = 1000,
#'     seed = 123456789,
#'     output_dir = "",
#'     output = "")
#' 
#' @author Vladislav BALAGURA <balagura@cern.ch>
#' @export
# [[Rcpp::export]]
sim <- function(n_points = 5000, n_turns = c(1000, 1000, 0, 5000),
                kick_model = "precise",
                n_sigma_cut = 5,
                density_and_field_interpolators_n_cells_along_grid_side = c(500, 500),
                n_random_points_to_check_interpolation = 10000,
                select_one_turn_out_of = 1000,
                seed = 123456789,
                output_dir = "tmp",
                output = "integrals.per.turn avr.xy.per.turn integrals.per.particle avr.xy.per.particle points") {
    list(n_points = n_points,
         n_turns = n_turns,
         kick_model = kick_model,
         n_sigma_cut = n_sigma_cut,
         density_and_field_interpolators_n_cells_along_grid_side =
             density_and_field_interpolators_n_cells_along_grid_side,
         n_random_points_to_check_interpolation = n_random_points_to_check_interpolation,
         select_one_turn_out_of = select_one_turn_out_of,
         seed = seed,
         output_dir = output_dir,
         output = output)
}
