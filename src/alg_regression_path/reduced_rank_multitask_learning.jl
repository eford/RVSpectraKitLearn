#=  ADMM Algorithmic Regularization Paths for Sparse Statistical Machine Learning
#   Yue Hu, Eric C. Chi and Genevera I. Allen
#   arXiv:1504.06637v1
# 3.2 Reduced-Rank Multi-Task Learning
  Goal:  fit for B in Y = X*B+eps
  n # num iid samples
  p # num covariates
  q # num outcomes
# Algorithm 4 Algorithmic Regularization Path for Reduced-Rank Regression
  1. Initialize: Z0 = 0, U0 = 0, g0 = e, and step size t > 0.
  2. Precompute: H = (XTX=n+I)􀀀1 and HXTY.
  3. While kZkk 6= 0
     gk = gk􀀀1 +t (or gk = gk􀀀1t).
     Bk = HXTY+H(Zk􀀀1􀀀Uk􀀀1).
     Zk = SVTgk (Bk+Uk􀀀1). (Record Zk at each iteration.)
  Uk = Uk􀀀1+Bk􀀀Uk
  end
  4. Output fZk : k = 1;    ;Kg as the algorithmic regularization path.
=#

# Singular Value Treshold of A w/ constant gammma
function SVT(A::Array,gamma::Float64)
    fac = svdfact(A)
    fac[:U]*(diagm(max(0,fac[:S].-gamma)))*fac[:Vt]
end

function fit_arp_rrr{T<:Number}(X::Array{T,2}, Y::Array{T,2}; stepsize::Float64 =1.05, gamma_init::Float64 =0.01, max_itterations::Integer=200, min_active::Integer=0, verbose::Bool=false)
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
  gamma = gamma_init
  last_itteration = max_itterations
  for k in 1:max_itterations
    gamma = gamma * stepsize
    B = HXtY+H*(Z-U)
    Z = SVT(B+U,gamma)
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
  # resize to last_itteration?
  return (gamma_list,Z_list)
end

#= Test orig alg
 n = 25
 p = 4
 q = 3
 X = randn(n,p)
 B = make_B(p,q,frac_active = 0.1)
 Y = X*B + 0.01*randn(n,q)
 println("B norms: ", length(find(B)), "  ", sumabs(B), " ", sumabs2(B))

(gamma_list, Z_list)  = fit_arp_rrr(X,Y, max_itterations = 200)
 print_path(gamma_list, Z_list)
=#
 