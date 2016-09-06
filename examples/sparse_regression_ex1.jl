using RvSpectraKitLearn

function make_B(p::Integer, q::Integer; frac_active::Float64 = 0.1, bias::Float64 = 0.0, sigma::Float64 = 1.0)
  B = zeros(p,q)
  for i in 1:p, j in 1:q
      if rand() < frac_active
        B[i,j] = bias+sigma*randn()
      end
  end
  return B
end

# Test sparse regression via  alg
n = 25
 p = 100
 q = 100
 frac_active = 0.1
 sigma_x = 0.01
 sigma_y = 0.01
 X = randn(n,p)
 B = make_B(p,q,frac_active=frac_active,bias=1.0,sigma=sigma_x)
 Y = X*B + sigma_y*randn(n,q)
 X_crossvalid = randn(n,p)
 Y_crossvalid = X_crossvalid*B + sigma_y*randn(n,q)
 println("B norms: ", length(find(B)), "  ", sumabs(B), " ", sumabs2(B))

(gamma_list, Z_list)  = fit_arp_slr_2d(X,Y,stepsize=1.05,max_itterations=1000, min_active=floor(Int64,p*q*frac_active*0.5))
(L0,L1,L2,rmsdelta) = compute_alr_path_stats(Z_list,target=B)
(chisq,chisq_refit,chisq_cv) = compute_alr_path_chisq!(Z_list,X,Y,X_crossvalid,Y_crossvalid)
print_alr_path(gamma_list, Z_list,X,Y)


#=
#Pkg.add("PyPlot")  # Run once if you don't already have PyPlot installed
using PyPlot

plt[:loglog](gamma_list,L0,"r-")
 plt[:loglog](gamma_list,L1/length(B),"b-")
 plt[:loglog](gamma_list,sqrt(L2),"g-")
 plt[:loglog](gamma_list,chisq_cv,"m-")
 plt[:loglog](gamma_list,chisq_cv,"m.")

plot(L0,chisq,"r-")
 plot(L0,chisq,"r.")
 plot(L0,chisq_refit,"b-")
 plot(L0,chisq_refit,"b.")
 plot(L0,chisq_cv,"g-")
 plot(L0,chisq_cv,"g.")

plot(L1,chisq,"r-")
 plot(L1,chisq,"r.")
 plot(L1,chisq_refit,"b-")
 plot(L1,chisq_refit,"b.")
 plot(L1,chisq_cv,"g-")
 plot(L1,chisq_cv,"g.")

plot(L0,L1./L0,"g-")
 plot(L0,L1./L0,"g.")

plot(L0,chisq,"r-")
 plot(L0,chisq,"r.")
 plot(L0,chisq_cv,"b-")
 plot(L0,chisq_cv,"b.")

=#
