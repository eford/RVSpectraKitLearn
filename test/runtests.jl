using RvSpectraKitLearn
using Base.Test

# write your own tests here
@test 1 == 1

include(joinpath(Pkg.dir("RvSpectraKitLearn"),"path_to_spectra.jl"))
cd(path_to_spectra)

(datafiles, phases) = get_filenames();
(lambda, obs) = read_filelist(datafiles);

# Test 1
using MultivariateStats
@test test_pca_eford(obs) < 1e-8

include("test2.jl")
