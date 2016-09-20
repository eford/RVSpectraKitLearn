using RvSpectraKitLearn
using MultivariateStats

# Likely need to update this to point to directory containing the spectra you wish to analyze
#include(joinpath(Pkg.dir("RvSpectraKitLearn"),"path_to_spectra.jl"))
path_to_spectra = joinpath(Pkg.dir(),"..","..","Code","XavierSpotData")
cd(path_to_spectra)

function get_filenames_multi_lats(; num_max_files::Integer = 0)
  files = readdir()
  datafiles = String[]
  strengths = Float64[]
  inclinations = Float64[]
  latitudes = Float64[]
  resolutions = Float64[]
  samplings = Float64[]
  phases = Float64[]
  for f in files
    m = match(r"(\w+)_(\d+)pc_incl(\d+)_lat(\d+)_reso(\d+)k_(-?\d+\.?\d*)pxsampl_(-?\d+\.?\d*)\.csv",f)
    if m==nothing continue end
    effect_name = m.captures[1]
    resolution = parse(m.captures[5])
    if resolution> 149
    push!(strengths,parse(m.captures[2]))
	  push!(inclinations,parse(m.captures[3]))
    push!(latitudes,parse(m.captures[4]))
    push!(resolutions,parse(m.captures[5]))
    push!(samplings,parse(m.captures[6]))
    push!(phases,parse(m.captures[7]))
    push!(datafiles,f)
    end
    if num_max_files>0 && length(datafiles)==num_max_files break end
  end
  @assert(length(datafiles)>0)
  return (datafiles,phases)
  #idx = sortperm(phases)
  #return (datafiles[idx],phases[idx])
end



(datafiles, phases) = get_filenames_multi_lats(); # num_max_files=4);
(lambda, obs_noisefree) = read_filelist(datafiles);

datafiles

sampling = 3.0
planet_amplitude = 10.0
rvs_true = 0.0*cos(2pi*phases)

snr = 200.0
obs = make_noisy_spectra(obs_noisefree,snr,sampling=sampling)

pca_out = fit_pca_default(obs);
pca_scores = transform(pca_out,obs)   # compute scores using Julia's PCA routine for comparison's sake
pca_eford_out = fit_pca_eford(obs);   # compute first few scores using Eric's itterative PCA routine

doppler_comp_simple = calc_doppler_component_simple(lambda,obs);
genpca_simple_out = fit_gen_pca_rv(obs,doppler_comp_simple)

rvs_out_simple = est_rvs_from_pc(obs,genpca_simple_out[1],genpca_simple_out[2][:,1])

genpca_simple_out_1 = genpca_simple_out[3]
genpca_simple_out_2 = genpca_simple_out[3]
plot(rvs_out_simple,genpca_simple_out[3][:,2],"r.")
plot(rvs_out_simple,genpca_simple_out[3][:,3],"b.")
plot(rvs_out_simple,genpca_simple_out[3][:,4],"g.")
 xlabel("Apparent Velocity (m/s)")
 ylabel("Next PC Coefficients")


plot(rvs_out_simple,genpca_simple_out[3][:,2],"r.")
lat10_mask = find(x->ismatch(r"lat10",x),datafiles)
lat40_mask = find(x->ismatch(r"lat40",x),datafiles)
lat70_mask = find(x->ismatch(r"lat70",x),datafiles)
plot(rvs_out_simple[lat70_mask],genpca_simple_out[3][lat70_mask,2],"b.")
plot(rvs_out_simple[lat40_mask],genpca_simple_out[3][lat40_mask,2],"g.")
plot(rvs_out_simple[lat10_mask],genpca_simple_out[3][lat10_mask,2],"m.")
 xlabel("Apparent Velocity (m/s)")
 ylabel("Next PC Coefficient")



abs(genpca_simple_out_1.-genpca_simple_out_2)

doppler_comp_gp = calc_doppler_component_gp(lambda,obs);
genpca_gp_out = fit_gen_pca_rv(obs,doppler_comp_gp)
rvs_out_gp = est_rvs_from_pc(obs,genpca_gp_out[1],genpca_gp_out[2][:,1])
rms_simple = sqrt(sumabs2(rvs_out_simple.-rvs_true)/length(rvs_true))
rms_gp = sqrt(sumabs2(rvs_out_gp.-rvs_true)/length(rvs_true))
println("snr = ", snr, " rms = ", rms_simple, "  ", rms_gp)

plot(rvs_out_simple,genpca_simple_out[3][:,3],"r.")
 plot(rvs_out_gp,genpca_gp_out[3][:,3],"b.")


using PyPlot
plot(phases,rvs_true,"b-")
 plot(phases,rvs_out_simple,"g.")
 xlabel("Phase")
 ylabel("Apparent Velocity (m/s)")

plot(phases,rvs_out_gp,"g.")
 plot(phases,rvs_out_simple.-rvs_true,"r.")
 plot(phases,rvs_out_gp.-rvs_true,"g.")

