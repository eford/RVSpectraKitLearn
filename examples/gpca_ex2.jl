using RvSpectraKitLearn
using MultivariateStats

# Likely need to update this to point to directory containing the spectra you wish to analyze
include(joinpath(Pkg.dir("RvSpectraKitLearn"),"path_to_spectra.jl"))
cd(path_to_spectra)

(datafiles, phases) = get_filenames(); # num_max_files=4);
(lambda, obs_noisefree) = read_filelist(datafiles);

sampling = 3.0
planet_amplitude = 10.0
rvs_true = -planet_amplitude*cos.(2pi*phases)

snr_list = [ 300.0, 200.0, 150.0, 100.0, 80.0, 60.0, 40.0 ]
rms_simple_list = Array(Float64,0)
rms_gp_list = Array(Float64,0)
for snr in snr_list
 obs = make_noisy_spectra(obs_noisefree,snr,sampling=sampling)
 # First, use simple finite differences on averaged spectra to compute derivatives and resulting RVs
 doppler_comp_simple = calc_doppler_component_simple(lambda,obs);
 genpca_simple_out = fit_gen_pca_rv(obs,doppler_comp_simple)
 rvs_out_simple = est_rvs_from_pc(obs,genpca_simple_out[1],genpca_simple_out[2][:,1])
 # Now, use GP to compute derivatives and resulting RVs
 doppler_comp_gp = calc_doppler_component_gp(lambda,obs);
 genpca_gp_out = fit_gen_pca_rv(obs,doppler_comp_gp)
 # Evaluate quality of measured RVs
 rvs_out_gp = est_rvs_from_pc(obs,genpca_gp_out[1],genpca_gp_out[2][:,1])
 rms_simple = sqrt(sumabs2(rvs_out_simple.-rvs_true)/length(rvs_true))
 rms_gp = sqrt(sumabs2(rvs_out_gp.-rvs_true)/length(rvs_true))
 push!(rms_simple_list,rms_simple)
 push!(rms_gp_list,rms_gp)
 println("snr = ", snr, " rms = ", rms_simple, "  ", rms_gp)
end

using PyPlot
plt[:loglog](snr_list,rms_simple_list,"r.")
 plt[:loglog](snr_list,rms_gp_list,"g.")
 xlabel("SNR per Resolution Element")
 ylabel("RMS Velocity (m/s)")
