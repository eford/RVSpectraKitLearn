#=  ADMM Algorithmic Regularization Paths for Sparse Statistical Machine Learning
#   Yue Hu, Eric C. Chi and Genevera I. Allen
#   arXiv:1504.06637v1
# 3.1 Sparse Linear Regression
  Goal:  fit for B in Y = X*B+eps
=#

function soft_threshold(x::Number, gamma::Float64)
  # Soft threshold of x to be writen to z
  gamma_by_2 = gamma/2
    if x > gamma_by_2
      return x - gamma_by_2
    elseif x < -gamma_by_2
      return x + gamma_by_2
    else
      return 0.0
    end
end

function soft_threshold!(x::Matrix, gamma::Float64, z::Matrix)
  @assert size(x) == size(z)
  for i in 1:size(x,1)
    for j in 1:size(x,2)
      z[i,j] = soft_threshold(x[i,j],gamma)
      #println("i=",i," j=",j," x=",x[i,j]," z=",z[i,j], " gamma=",gamma)
    end
  end
  return z
end

function fit_arp_slr_2d{T<:Number}(X::Array{T,2}, Y::Array{T,2}; stepsize::Float64 =1.05, gamma_init::Float64 =0.01, max_itterations::Integer=200, min_active::Integer=0, verbose::Bool=false)
  @assert size(X,1) == size(Y,1)
  n = size(X,1)
  p = size(X,2)
  q = size(Y,2)
  gamma_list = zeros(T,max_itterations)
  chisq_list = zeros(T,max_itterations)
  Z_list = zeros(max_itterations,p,q)
  B = zeros(p,q)  # Preallocate working matrices
  Z = zeros(p,q)
  U = zeros(p,q)
  H = inv(X'*X/n+diagm(ones(p)))
  HXtY = H*X'*Y
  local gamma = gamma_init
  last_itteration = max_itterations
  #println("# fit_arp_slr_2d")
  for k in 1:max_itterations
    gamma = gamma * stepsize
    B = HXtY+H*(Z-U)
    soft_threshold!(B+U,gamma,Z)
    gamma_list[k]= gamma
    Z_list[k,:,:] = Z
    U = U + B - Z
    num_active = length(find(Z))
    if verbose
      chisq = sumabs2(Y-X*Z)
      print("# ", k, " ", gamma,": ")
      print("Z norms= ", num_active, " ", sumabs(Z), " ", sumabs2(Z), " chi^2= ", chisq, " ")
      #print("B=",B," Z=",Z," U=",U," ")
      println("")
    end
    if sumabs2(Z) == 0. || num_active < min_active
      last_itteration = k
      break
    end
  end
  resize!(gamma_list,last_itteration)
  resize!(chisq_list,last_itteration)
  # resize to last_itteration?
  return (gamma_list,Z_list[1:last_itteration,:,:])
end

# Simple linear regression to update the non-zero coefficients
function refit_arp_slr_2d{T<:Number}(X::Array{T,2}, Y::Array{T,2}, B::Array{T,2})
  @assert size(X,1) == size(Y,1)
  p = size(X,2)
  Z = copy(B)
  for col in 1:size(B,2)         # Separate mutliple parameter regression problems into 1 at a time
    idx_active = find(B[:,col])  # find non-zero entries
    Xact = X[:,idx_active]
    Bact = (Xact'*Xact) \ (Xact'*Y[:,col])
    Z[idx_active,col] = vec(Bact)
  end
  pred = X*Z
  chisq = sumabs2(Y-pred)
  return ( Z,chisq )
end

function refit_arp_slr_2d{T<:Number}(X::Array{T,2}, Y::Array{T,2}, B::Array{T,2}, Xcv::Array{T,2}, Ycv::Array{T,2} )
  @assert size(X,1) == size(Y,1)
  @assert size(X) == size(Xcv)
  @assert size(Y) == size(Ycv)
  (Z, chisq_refit) = refit_arp_slr_2d(X,Y,B)
  chisq_cv = sumabs2(Ycv-Xcv*Z)
  return (Z,chisq_refit, chisq_cv )
end


#=
using Optim

function refit_arp_slr_2d_nonlin{T<:Number}(X::Array{T,2}, Y::Array{T,2}, B::Array{T,2})
  @assert size(X,1) == size(Y,1)
  p = size(X,2)

  # Function to minimize
  score(Z::Array{T,2}) = sumabs2(Y-X*Z)
  idx_active = find(B)
  function score(param::Array{Float64,1})
     Z = copy(B)
     Z[idx_active] = param
     score(Z)
  end
  # This is much slower than linear regression.  Left code in case the problem becomes weakly non-linear in future
  opt_result = optimize(score,B[idx_active],LBFGS())
  Z = copy(B)
  Z[idx_active] = Optim.minimizer(opt_result)
  return ( Z,Optim.minimum(opt_result) )
end
=#

