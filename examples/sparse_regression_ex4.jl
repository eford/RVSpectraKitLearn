using RvSpectraKitLearn

# Likely need to update this to point to directory containing the spectra you wish to analyze
include(joinpath(Pkg.dir("RvSpectraKitLearn"),"path_to_spectra.jl"))
cd(path_to_spectra)
(datafiles_planet, phases_planet) = get_filenames();
(lambda_planet, obs_planet_noisefree) = read_filelist(datafiles_planet);
planet_amplitude = 10.0
rvs_planet = -planet_amplitude*cos(phases_planet*2pi)

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

(datafiles_spot, phases_spot) = get_filenames_multi_lats(); # num_max_files=4);
(lambda_spot, obs_spot_noisefree) = read_filelist(datafiles_spot);
rvs_spot = zeros(phases_spot)

rvs_all = [rvs_spot, rvs_planet]
@assert(all(lambda_spot .== lambda_planet))
obs_noisefree = Array(Float64,(size(obs_planet_noisefree,1),size(obs_planet_noisefree,2)+size(obs_spot_noisefree,2)))
obs_noisefree[:,1:size(obs_spot_noisefree,2)] = obs_spot_noisefree
obs_noisefree[:,1+size(obs_spot_noisefree,2):end] = obs_planet_noisefree
obs_noisefree

function find_important_lamba{T}(lambda::Array{T,1}, obs::Array{T,2}, rvs::Array{T,1};
                                 chunk_size::Integer=1500, frac_to_keep = 0.1 )
  Y = reshape(rvs,(length(rvs),1))
 beta = ones(length(lambda)) # coefficients
 new_beta = copy(beta)
  idx_active = find(beta)
  num_chunks = ceil(Int64,length(idx_active)/chunk_size)
  for c in 1:num_chunks
    idx_min = 1+(c-1)*chunk_size
    idx_max = min(c*chunk_size,length(idx_active))
    lambda_range = idx_active[idx_min:idx_max]
    #println("idx=",idx_min,"-",idx_max,"  active=",idx_active[idx_min],"-",idx_active[idx_max]," lambda=",lambda[idx_active[idx_min]],"-",lambda[idx_active[idx_max]])
    X = obs[lambda_range,:]'
    (gamma_list, Z_list)  = fit_arp_slr_2d(X,Y,stepsize=1.01,gamma_init=0.001, max_itterations=1000, min_active=floor(Int64,size(X,2)*size(Y,2)*frac_to_keep) )
    new_beta[lambda_range] = vec(Z_list[size(Z_list,1),:,1])
    #println("size(Z_list)=", size(Z_list,1), " nactive=",length(find(new_beta[lambda_range])))
  end
  beta[:] = new_beta
 return new_beta
end

beta = find_important_lamba(lambda_spot,obs_noisefree,rvs_all, chunk_size=200, frac_to_keep = [0.1] )

#Pkg.add("JLD")  # If you don't already have JLD installed
using JLD
save("important_lambda_0.1.jld", "frac_to_keep", [0.1], "beta", beta)

#=
# Plot which wavelengths are being used from the last itteration of the algorithm regularization path algorithm
=#
using PyPlot
lambda = lambda_planet
obs = obs_noisefree
lambda_range = 1:length(lambda_planet)
 zmax = maximum(abs(beta[lambda_range]))
  ymax = maximum(vec(obs[lambda_range,1]))
  ymin = minimum(vec(obs[lambda_range,1]))
  yrange = ymax-ymin
  #plot(lambda,vec(obs[:,1]),"r.")
  min_lambda = lambda[lambda_range][1]
  max_lambda = lambda[lambda_range][end]
  idx_around = find(x->min_lambda<=x<=max_lambda,lambda)
  plot(lambda[idx_around],vec(obs[idx_around,1]),"r-")
  plot(lambda[lambda_range],ymin+yrange/2*(1+1/zmax*(beta[lambda_range])),"b+")
  plot(lambda[lambda_range],ymin+yrange/2*(ones(length(lambda_range))),"b-")

