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

# Test eford alg
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
 #print_path(gamma_list, Z_list,X,Y)
 (L0,L1,L2,rmsdelta) = compute_alr_path_stats(Z_list,target=B)
 (chisq,chisq_refit,chisq_cv) = compute_alr_path_chisq(Z_list,X,Y,X_crossvalid,Y_crossvalid)
#(chisq,chisq_refit,chisq_cv) = compute_alr_path_chisq(Z_list[1:25,:,:],X,Y,X_crossvalid,Y_crossvalid)
#chisq_cv = compute_alr_path_chisq(Z_list,X_crossvalid,Y_crossvalid)

keep_path = find(x->0.5*length(find(B))<=x<=1.3*length(find(B)),L0)
Z_list[keep_path,:,:]

for i in 1:4:length(keep_path)
  j = keep_path[i]
  Z = reshape(Z_list[j,:,:],(size(Z_list,2),size(Z_list,3)))
  chisq_alr = chisq[j] /(length(Y)*sigma_y^2)
  (Zr,rchisq_refit,rchisq_cv) = refit_arp_slr_2d(X,Y,Z,X_crossvalid,Y_crossvalid)
  rchisq_refit /= (length(Y)*sigma_y^2)
  rchisq_cv  /= (length(Y)*sigma_y^2)
  @printf "%d %5.2f: %4d %8g  X^2= %8g %8g %8g\n" j gamma_list[j] L0[j] L1[j] rchisq_alr  rchisq_refit  rchisq_cv
end

using PyPlot

plot(L0,chisq,"r-")
 plot(L0,chisq,"r.")
 plot(L0,chisq_refit,"b-")
 plot(L0,chisq_refit,"b.")
 plot(L0,chisq_cv,"g-")
 plot(L0,chisq_cv,"g.")
 #plt[:loglog](L0,gamma_list,"g.")

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

plt[:loglog](gamma_list,L0,"r-")
 plt[:loglog](gamma_list,L1/length(B),"b-")
 plt[:loglog](gamma_list,sqrt(L2),"g-")
 plt[:loglog](gamma_list,chisq_cv,"m-")
 plt[:loglog](gamma_list,chisq_cv,"m.")

 plt[:loglog](gamma_list,sqrt(sigma_x*L2/length(B)),"g-")
 plt[:loglog](gamma_list,sqrt(chisq*sigma_y^2/(length(Y))),"r+")
 plt[:loglog](gamma_list,sqrt(chisq_cv*sigma_y^2/(length(Y))),"b+")
 plt[:loglog](gamma_list,rmsdelta,"g+")

