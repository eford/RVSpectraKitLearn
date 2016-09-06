using RvSpectraKitLearn

# Likely need to update this to point to directory containing the spectra you wish to analyze
include(joinpath(Pkg.dir("RvSpectraKitLearn"),"path_to_spectra.jl"))
cd(path_to_spectra)

(datafiles, phases) = get_filenames();
(lambda, obs) = read_filelist(datafiles);

lambda_min = 4500
 lambda_max = 4510
 idx_min = findfirst(x->lambda_min<x<lambda_max,lambda)
 idx_max = findlast(x->lambda_min<x<lambda_max,lambda)
 lambda_range = idx_min:idx_max

Y = reshape(-10.0*cos(phases*2pi),(length(phases),1))
 X = obs[lambda_range,:]'
 @assert size(X,1) == size(Y,1)
 n = size(X,1)
 p = size(X,2)
 q = size(Y,2)

@time (gamma_list, Z_list)  = fit_arp_slr_2d(X,Y,stepsize=1.005,gamma_init=0.001, max_itterations=1000, min_active=floor(Int64,p*q/3))
 (L0,L1,L2,rmsdelta) = compute_alr_path_stats(Z_list)
 (chisq,chisq_refit) = compute_alr_path_chisq(Z_list,X,Y)
 print_alr_path(gamma_list, Z_list,X,Y)
 Zr_list = similar(Z_list)


# Check what wavelengths are being used from the last itteration of the algorithm regularization path algorithm
#using PyPlot
Zr_idx = size(Z_list,1)
 Z = reshape(Z_list[Zr_idx,:,:],(size(Z_list,2),size(Z_list,3)))
  rchisq_alr = chisq[Zr_idx] /(length(Y))
  @time (Zr,rchisq_refit) = refit_arp_slr_2d(X,Y,Z)
  Zr_list[Zr_idx,:,:] = Zr
  pred = X*reshape(Zr,(p,q))
  find(reshape(Zr,(p,q)))
  plot(lambda[lambda_range],vec(X[1,:]),"r.")
  plot(lambda[lambda_range],3000+10000*(reshape(Zr,(p,q))),"b+")

idx_active = find(Z[:,1])
Xact = X[:,idx_active]
X[:,idx_active]*Z[idx_active,1]
Bact = (Xact'*Xact) \ (Xact'*Y)
Y-X[:,idx_active]*Bact


# Check that the predictions are good
plot(phases,Y,"r.")
 plot(phases,Y,"r+")
 plot(phases,pred,"b+")


# Refit for several values of regularization parameter
Zr_list = similar(Z_list)
for i in 1:10:size(Z_list,1)
  Z = reshape(Z_list[i,:,:],(size(Z_list,2),size(Z_list,3)))
  rchisq_alr = chisq[i] /(length(Y))
  (Zr,rchisq_refit) = refit_arp_slr_2d(X,Y,Z)
  Zr_list[i,:,:] = Zr
  rchisq_refit /= (length(Y))
  @printf "%d %5.2f: %4d %8g  X^2= %8g %8g\n" i gamma_list[i] L0[i] L1[i] rchisq_alr  rchisq_refit
end


# Random other diagnostic plots
plot(L0,chisq,"r-")
 plot(L0,chisq,"r.")
 plot(L0,chisq_refit,"b-")
 plot(L0,chisq_refit,"b.")
 #plot(L0,chisq_cv,"g-")
 #plot(L0,chisq_cv,"g.")

plt[:loglog](L0,gamma_list,"g.")

plot(L1,chisq,"r-")
 plot(L1,chisq,"r.")
 plot(L1,chisq_refit,"b-")
 plot(L1,chisq_refit,"b.")

plot(L0,L1./L0,"g-")
 plot(L0,L1./L0,"g.")

plot(L0,chisq,"r-")
 plot(L0,chisq,"r.")

plt[:loglog](gamma_list,L0,"r-")
 plt[:loglog](gamma_list,L1/(p*q),"b-")
 plt[:loglog](gamma_list,sqrt(L2),"g-")

 plt[:loglog](gamma_list,sqrt(L2/(p*q)),"g+")
 plt[:loglog](gamma_list,sqrt(chisq/(n*q)),"r+")
 plt[:loglog](gamma_list,rmsdelta,"g+")

=#
