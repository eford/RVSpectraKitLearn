using MultivariateStats

# Simply calling Julia's PCA routine for testing purposes
function fit_pca_default(X::Array; max_num_components::Integer=5)
  num_lambda = size(X,1)
  num_spectra = size(X,2)
  M = fit(PCA,X,maxoutdim=max_num_components,pratio=0.999999)
  num_components = size(projection(M),2)
  #Y = transform(M,X)
  #Xr = reconstruct(M,Y)
  #save("pca.jld", "num_lambda", num_lambda, "num_spectra", num_spectra, "num_components", num_components, "X", X, "M", M, "Y", Y, "Xr", Xr)
  return M
end

# Compute the PCA component with the largest eigenvalue
# X is data, r is vector of random numbers, s is preallocated memory; r & s  are of same length as each data point
function compute_pca_component!{T}(X::AbstractArray{T,2}, r::AbstractArray{T,1}, s::AbstractArray{T,1}; tol::Float64=1e-8, max_it::Int64 = 10 )
  num_lambda = size(X,1)
  num_spectra = size(X,2)
  @assert(length(r)==num_lambda)
  #rand!(r)                                               # assume r is already randomized
  last_mag_s = 0.0
  #s = zeros(T,num_lambda)
  for c in 1:max_it
	s[:] = zero(T)
	for i in 1:num_spectra
	   BLAS.axpy!(dot(view(X,:,i),r),view(X,:,i),s)         # s += dot(X[:,i],r)*X[:,i]
	end
	mag_s = sqrt(sumabs2(s))
	r[:]  = s/mag_s
	if abs(mag_s-last_mag_s) < tol*mag_s break end
	last_mag_s = mag_s
  end
  return r
end

# Compute the PCA component with the largest eigenvalue
function compute_pca_component{T}(X::AbstractArray{T,2}; tol::Float64=1e-12, max_it::Int64 = 10)
  num_lambda = size(X,1)
  r = rand(T,num_lambda)
  s = zeros(T,num_lambda)
  compute_pca_component!(X,r,s,tol=tol,max_it=max_it)
end

# Compute first num_components basis vectors for PCA, only return the basis vectors
function fit_pca_eford_proj_only{T}(X::Array{T,2}; num_components::Integer=4, tol::Float64=1e-12, max_it::Int64 = 10)
  num_lambda = size(X,1)
  num_spectra = size(X,2)
  M = rand(T, (num_lambda,num_components) )  # random initialization is part of algorithm (i.e., not zeros)
  s = zeros(T,num_lambda)                    # pre-allocated memory for compute_pca_component
  Xtmp = X.-mean(X,2)                        # perform PCA after subtracting off mean
  compute_pca_component!(Xtmp,view(M,:,1),s,tol=tol)  # compute first component separately, since generalized algorithm will need to
  for j in 2:num_components
     for i in 1:num_spectra
       score = dot(view(Xtmp,:,i),view(M,:,j-1))
       Xtmp[:,i] .-= score*view(M,:,j-1)
	   end
	  compute_pca_component!(Xtmp,view(M,:,j),s,tol=tol,max_it=max_it)
  end
  return M
end

# Compute first num_components basis vectors for PCA, also return scores and print even more info
function fit_pca_eford{T}(X::Array{T,2}; num_components::Integer=4, tol::Float64=1e-12, max_it::Int64 = 10)
  num_lambda = size(X,1)
  num_spectra = size(X,2)
  M = rand(T, (num_lambda,num_components) )  # random initialization is part of algorithm (i.e., not zeros)
  s = zeros(T,num_lambda)                    # pre-allocated memory for compute_pca_component
  scores = zeros(num_spectra,num_components)
  mu = vec(mean(X,2))
  Xtmp = X.-mu                               # perform PCA after subtracting off mean
  totalvar = sumabs2(Xtmp)
  #println("# j = ", 0, " sumabs2(Xtmp) = ", totalvar)
  compute_pca_component!(Xtmp,view(M,:,1),s,tol=tol)  # compute first component separately, since generalized algorithm will need to
  fracvar = zeros(num_components)
  compute_pca_component!(Xtmp,view(M,:,1),s,tol=tol)
  for j in 2:num_components
	  for i in 1:num_spectra
      scores[i,j-1] = dot(view(Xtmp,:,i),view(M,:,j-1)) #/sumabs2(view(M,:,j-1))
      Xtmp[:,i] .-= scores[i,j-1]*view(M,:,j-1)
	  end
    fracvar[j-1] = sumabs2(Xtmp)/totalvar
    println("# j = ", j-1, " sumabs2(Xtmp) = ", sumabs2(Xtmp), " frac_var_remain= ", fracvar[j-1] )
	  compute_pca_component!(Xtmp,view(M,:,j),s,tol=tol,max_it=max_it)
  end
	for i in 1:num_spectra
       scores[i,num_components] = dot(view(Xtmp,:,i),view(M,:,num_components))
       Xtmp[:,i] .-= scores[i,num_components]*view(M,:,num_components)
	end
  fracvar[num_components] = sumabs2(Xtmp)/totalvar
  println("# j = ", num_components, " sumabs2(Xtmp) = ", sumabs2(Xtmp), " fracvar= ", fracvar[num_components] )
	return (mu,M,scores)
end

# test that eford's PCA routine is giving accurate results for dataset X and min(num_components,number of components kept by Julia's PCA routine)
function test_pca_eford{T}(obs::Array{T,2}; num_components::Integer=4 )
  pca_out = fit_pca_default(obs)
  (pca_eford_mu, pca_eford_out, pca_eford_scores) = fit_pca_eford(obs,num_components=num_components)
  for i in 1:min(size(projection(pca_out),2),size(pca_eford_out,2))
     if ! (projection(pca_out)[1,i] * pca_eford_out[1,i] > 0.0)
	    pca_eford_out[:,i] = -pca_eford_out[:,i]
        pca_eford_scores[:,i] = -pca_eford_scores[:,i]
	 end
  end
  maxabsdif = maximum(abs(projection(pca_out).-pca_eford_out[:,1:size(projection(pca_out),2)]))
  return maxabsdif
end



