module RvSpectraKitLearn

export fit_gen_pca_rv, est_rvs_from_pc
include("alt_pca.jl")

export calc_doppler_component_simple, calc_doppler_quadratic_term_simple
include("deriv_spectra_simple.jl")
#include("deriv_spectra_gp.jl")

export fit_pca_default, fit_pca_eford_proj_only, fit_pca_eford, test_pca_eford
include("generalized_pca.jl")

export get_filenames, read_filelist
include("read_inputs.jl")
#include("rand_junk.jl")

end # module



