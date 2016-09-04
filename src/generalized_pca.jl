
# Compute first num_components basis vectors for PCA, after subtracting projection onto fixed_comp
function fit_gen_pca_rv{T}(X::Array{T,2}, fixed_comp::Array{T,1}; quad_term::Array{T,1} = zeros(T,0), num_components::Integer=4, tol::Float64=1e-12, max_it::Int64 = 10)
  num_lambda = size(X,1)
  num_spectra = size(X,2)
  M = rand(T, (num_lambda,num_components) )  # random initialization is part of algorithm (i.e., not zeros)
  s = zeros(T,num_lambda)                    # pre-allocated memory for compute_pca_component
  scores = zeros(num_spectra,num_components)
  mu = vec(mean(X,2))
  Xtmp = X.-mu                               # perform PCA after subtracting off mean
  totalvar = sumabs2(Xtmp)
  fracvar = zeros(num_components)
  M[:,1] = fixed_comp                        # Force fixed (i.e., Doppler) component to replace first PCA component
  fixed_comp_norm = 1/sumabs2(fixed_comp)
  for i in 1:num_spectra
    scores[i,1] = z = (dot(view(Xtmp,:,i),fixed_comp)*fixed_comp_norm)  # Normalize differently, so scores are z (i.e., doppler shift)
	  Xtmp[:,i] -= z*fixed_comp
    if length(quad_term) == length(fixed_comp)  # Optionally subtract off quadratic term from Doppler shift
      Xtmp[:,i] -= 0.5*z*z*quad_term   # subtract off projection onto fixed component and quadratic term separately
    end
  end
  fracvar[1] = sumabs2(Xtmp)/totalvar
  #println("# sumabs2(Xtmp) = ", sumabs2(Xtmp) )
  #println("# j = ", 1, " sumabs2(Xtmp) = ", sumabs2(Xtmp), " frac_var_remain= ", fracvar[1] )
  for j in 2:num_components
    compute_pca_component!(Xtmp,view(M,:,j),s,tol=tol,max_it=max_it)
    for i in 1:num_spectra
       scores[i,j] = dot(view(Xtmp,:,i),view(M,:,j)) #/sumabs2(view(M,:,j-1))
       Xtmp[:,i] .-= scores[i,j]*view(M,:,j)
	  end
    fracvar[j] = sumabs2(Xtmp)/totalvar
    println("# j = ", j, " sumabs2(Xtmp) = ", sumabs2(Xtmp), " frac_var_remain= ", fracvar[j] )
  end
  return (mu,M,scores)
end

function est_rvs_from_pc{T}(X::Array{T,2}, mu::Vector{T}, doppler_comp::Vector{T})
  num_lambda = size(X,1)
  num_spectra = size(X,2)
  speed_of_light = 3.0e8
  rv = zeros(num_spectra)
  Xtmp = X.-mu
  fixed_comp_norm = 1/sumabs2(doppler_comp)
  for i in 1:num_spectra
      rv[i] = dot(Xtmp[:,i],doppler_comp)*fixed_comp_norm*speed_of_light
  end
  rv
end


