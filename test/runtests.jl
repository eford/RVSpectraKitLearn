using RvSpectraKitLearn
using Base.Test

# write your own tests here
@test 1 == 1

# Replace path for your system, or just comment out and run from relevant directory
path_to_spectra = "C:\\Users\\eford\\Box Sync\\SOAP simulations\\Aug2016_workshop\\SOAP_Spectra\\planet_10ms_150k"
cd(path_to_spectra)  

(datafiles, phases) = get_filenames();
(lambda, obs) = read_filelist(datafiles);

# Test 1
using MultivariateStats
@test test_pca_eford(obs) < 1e-8

include("test2.jl")