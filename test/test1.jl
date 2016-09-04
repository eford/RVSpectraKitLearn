using RvSpectraKitLearn

# Replace path for your system, or just comment out and run from relevant directory
path_to_spectra = "C:\\Users\\eford\\Box Sync\\SOAP simulations\\Aug2016_workshop\\SOAP_Spectra\\planet_10ms_150k"
cd(path_to_spectra)  

(datafiles, phases) = get_filenames();
(lambda, obs) = read_filelist(datafiles);
pca_out = fit_pca_default(obs);

using MultivariateStats
pca_scores = transform(pca_out,obs)   # compute scores using Julia's PCA routine for comparison's sake
pca_eford_out = fit_pca_eford(obs);
max_princ_comp_diff = test_pca_eford(obs)
@assert max_princ_comp_diff < 1e-8

doppler_comp_simple = calc_doppler_component_simple(lambda,obs);
genpca_simple_out = fit_gen_pca_rv(obs,doppler_comp_simple)
rvs_out = est_rvs_from_pc(obs,genpca_simple_out[1],genpca_simple_out[2][:,1])

# Untested, but shows outline for how would improve accuracy w/ quadratic term
#doppler_quadratic_term_simple = calc_doppler_quadratic_term_simple(lambda,obs);
#genpca_simple_quad_out = fit_gen_pca_rv(obs,doppler_comp_simple,quad_term=doppler_quadratic_term_simple)
#est_rvs_from_pc(obs,genpca_simple_out[1],genpca_simple_out[2][:,1])


