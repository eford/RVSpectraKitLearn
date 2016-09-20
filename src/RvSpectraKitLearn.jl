module RvSpectraKitLearn

using Compat
import Compat.view

export fit_gen_pca_rv, est_rvs_from_pc
include("alt_pca.jl")

export fit_pca_default, fit_pca_eford_proj_only, fit_pca_eford, test_pca_eford
include("generalized_pca.jl")

export calc_doppler_component_simple, calc_doppler_quadratic_term_simple
include("deriv_spectra_simple.jl")
export calc_doppler_component_gp
include("deriv_spectra_gp.jl")

export fit_arp_slr_2d, refit_arp_slr_2d
include("alg_regression_path/sparse_linear_regression.jl")
export print_alr_path, compute_alr_path_stats, compute_alr_path_chisq!
include("alg_regression_path/util.jl")

#export fit_arp_rrr
#include("alg_regression_path/reduced_rank_multitask_learning.jl")  # Implemented, but untested

export get_filenames, read_filelist, make_noisy_spectrum, make_noisy_spectra
include("read_inputs.jl")
#include("rand_junk.jl")

end # module



