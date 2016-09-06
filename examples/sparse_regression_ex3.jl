using RvSpectraKitLearn

# Likely need to update this to point to directory containing the spectra you wish to analyze
include(joinpath(Pkg.dir("RvSpectraKitLearn"),"path_to_spectra.jl"))
cd(path_to_spectra)

(datafiles, phases) = get_filenames();
(lambda, obs) = read_filelist(datafiles);

function find_important_lamba{T}(lambda::Array{T,1}, obs::Array{T,2}, phases::Array{T,1};
                                 chunk_size::Integer=1500, frac_to_keep = [0.2,0.5] )
  Y = reshape(-10.0*cos(phases*2pi),(length(phases),1))
 beta = ones(length(lambda)) # coefficients
 new_beta = copy(beta)
 for i in 1:length(frac_to_keep)
  idx_active = find(beta)
  println("i=",i," len=",length(idx_active))
  num_chunks = ceil(Int64,length(idx_active)/chunck_size)
  for c in 1:num_chunks
    idx_min = 1+(c-1)*chunck_size
    idx_max = min(c*chunck_size,length(idx_active))
    lambda_range = idx_active[idx_min:idx_max]
    #println("idx=",idx_min,"-",idx_max,"  active=",idx_active[idx_min],"-",idx_active[idx_max]," lambda=",lambda[idx_active[idx_min]],"-",lambda[idx_active[idx_max]])
    X = obs[lambda_range,:]'
    (gamma_list, Z_list)  = fit_arp_slr_2d(X,Y,stepsize=1.01,gamma_init=0.001, max_itterations=1000, min_active=floor(Int64,size(X,2)*size(Y,2)*frac_to_keep[i]) )
    new_beta[lambda_range] = vec(Z_list[size(Z_list,1),:,1])
    #println("size(Z_list)=", size(Z_list,1), " nactive=",length(find(new_beta[lambda_range])))
  end
  beta[:] = new_beta
 end
 return (lambda_range, beta)
end

(lambda_range, beta) = find_important_lamba(lambda,obs,phases,chunk_size=1500,frac_to_keep=[0.25,0.4])
#=
# Plot which wavelengths are being used from the last itteration of the algorithm regularization path algorithm
=#
using PyPlot

zmax = maximum(abs(beta))
  ymax = maximum(vec(obs[:,1]))
  ymin = minimum(vec(obs[:,1]))
  yrange = ymax-ymin
  #plot(lambda,vec(obs[:,1]),"r.")
  min_lambda = lambda[lambda_range][1]
  max_lambda = lambda[lambda_range][end]
  idx_around = find(x->min_lambda<=x<=max_lambda,lambda)
  plot(lambda[idx_around],vec(obs[idx_around,1]),"r-")
  plot(lambda[lambda_range],ymin+yrange/2*(1+1/zmax*(beta[lambda_range])),"b+")
  plot(lambda[lambda_range],ymin+yrange/2*(ones(length(lambda_range))),"b-")

#=
=#

#=
# Check that the predictions are good
pred = X*reshape(Z,(p,q))
 plot(phases,Y,"r.")
 plot(phases,Y,"r+")
 plot(phases,pred,"b+")
=#

#=
# Random other diagnostic plots
plot(L0,chisq,"r-")
 plot(L0,chisq,"r.")
 plot(L0,chisq_refit,"b-")
 plot(L0,chisq_refit,"b.")
 #plot(L0,chisq_cv,"g-")
 #plot(L0,chisq_cv,"g.")

plt[:loglog](L0,gamma_list,"g.")
 plot(L0,L1./L0,"g-")
 plot(L0,L1./L0,"g.")

plot(L1,chisq,"r-")
 plot(L1,chisq,"r.")
 plot(L1,chisq_refit,"b-")
 plot(L1,chisq_refit,"b.")


plot(L0,chisq,"r-")
 plot(L0,chisq,"r.")

plt[:loglog](gamma_list,L0,"r-")
 plt[:loglog](gamma_list,L1/(p*q),"b-")
 plt[:loglog](gamma_list,sqrt(L2),"g-")

 plt[:loglog](gamma_list,sqrt(L2/(p*q)),"g+")
 plt[:loglog](gamma_list,sqrt(chisq/(n*q)),"r+")
 plt[:loglog](gamma_list,rmsdelta,"g+")

=#
