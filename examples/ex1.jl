using RvSpectraKitLearn
using MultivariateStats

# Likely need to update this to point to directory containing the spectra you wish to analyze
include(joinpath(Pkg.dir("RvSpectraKitLearn"),"path_to_spectra.jl"))
cd(path_to_spectra)  

(datafiles, phases) = get_filenames();
(lambda, obs) = read_filelist(datafiles);

pca_out = fit_pca_default(obs);
pca_scores = transform(pca_out,obs)   # compute scores using Julia's PCA routine for comparison's sake
pca_eford_out = fit_pca_eford(obs);   # compute first few scores using Eric's itterative PCA routine

doppler_comp_simple = calc_doppler_component_simple(lambda,obs);
genpca_simple_out = fit_gen_pca_rv(obs,doppler_comp_simple)
rvs_out = est_rvs_from_pc(obs,genpca_simple_out[1],genpca_simple_out[2][:,1])
