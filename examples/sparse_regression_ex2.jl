using RvSpectraKitLearn

# Likely need to update this to point to directory containing the spectra you wish to analyze
include(joinpath(Pkg.dir("RvSpectraKitLearn"),"path_to_spectra.jl"))
cd(path_to_spectra)

(datafiles, phases) = get_filenames();
(lambda, obs) = read_filelist(datafiles);

lambda_min = 4500
 lambda_max = 4520
 idx_min = findfirst(x->lambda_min<x<lambda_max,lambda)
 idx_max = findlast(x->lambda_min<x<lambda_max,lambda)
 lambda_range = idx_min:idx_max

Y = reshape(-10.0*cos(phases*2pi),(length(phases),1))
 X = obs[lambda_range,:]'
 @assert size(X,1) == size(Y,1)
 n = size(X,1)
 p = size(X,2)
 q = size(Y,2)

@time (gamma_list, Z_list)  = fit_arp_slr_2d(X,Y,stepsize=1.01,gamma_init=0.001, max_itterations=1000, min_active=floor(Int64,p*q/10))
(L0,L1,L2,rmsdelta) = compute_alr_path_stats(Z_list)
#(chisq,chisq_refit) = compute_alr_path_chisq!(Z_list,X,Y)  # Refitting via least squares makes coefficients worse, since highly underdetermined
print_alr_path(gamma_list, Z_list,X,Y)

#=
# Plot which wavelengths are being used from the last itteration of the algorithm regularization path algorithm
=#
using PyPlot
Z_idx = size(Z_list,1) # floor(Int64,size(Z_list,1)/2)
  Z = reshape(Z_list[Z_idx,:,:],(size(Z_list,2),size(Z_list,3)))
  find(reshape(Z,(p,q)))
  zmax = maximum(Z)
  ymax = maximum(vec(X[1,:]))
  ymin = minimum(vec(X[1,:]))
  yrange = ymax-ymin
  plot(lambda[lambda_range],vec(X[1,:]),"r.")
  plot(lambda[lambda_range],ymin+yrange/2*(1+1/zmax*(reshape(Z,(p,q)))),"b+")
  plot(lambda[lambda_range],ymin+yrange/2*ones(length(lambda_range)),"-")
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
