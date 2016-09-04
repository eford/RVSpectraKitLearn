
doppler_comp_simple = calc_doppler_component_simple(lambda,obs);
genpca_simple_out = fit_gen_pca_rv(obs,doppler_comp_simple)
rvs_out = est_rvs_from_pc(obs,genpca_simple_out[1],genpca_simple_out[2][:,1])

true_rvs = [-10.0*cos(2pi*phases[i]) for i in 1:25]  # not actually sure what true rvs are 

@test maximum(abs(rvs_out .- true_rvs)) < 0.1
