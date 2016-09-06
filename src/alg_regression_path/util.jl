function make_B(p::Integer, q::Integer; frac_active::Float64 = 0.1, bias::Float64 = 0.0, sigma::Float64 = 1.0)
  B = zeros(p,q)
  for i in 1:p, j in 1:q
      if rand() < frac_active
        B[i,j] = bias+sigma*randn()
      end
  end
  return B
end

function print_alr_path{T<:Number}(gamma_list::Vector,Z_list::Array{T,3})
  @assert length(gamma_list) == size(Z_list,1)
  for i in 1:length(gamma_list)
    Z = Z_list[i,:,:]
    num_zeros = length(Z)-length(find(Z))
    L0_Z = num_nonzeros = length(find(Z))
    L1_Z = sumabs(Z)
    L2_Z = sumabs2(Z)
    #rank_Z = rank(Z)
    println(i, "  ", gamma_list[i], ": ", L0_Z, "  ", L1_Z, " ", L2_Z)# " rank= ", rank_Z, " Z=",Z)
  end
end

function print_alr_path{T<:Number}(gamma_list::Vector,Z_list::Array{T,3},X::Array{T,2},Y::Array{T,2})
  @assert length(gamma_list) == size(Z_list,1)
  for i in 1:length(gamma_list)
    Z = reshape(Z_list[i,:,:],(size(Z_list,2),size(Z_list,3)))
    num_zeros = length(Z)-length(find(Z))
    L0_Z = num_nonzeros = length(find(Z))
    L1_Z = sumabs(Z)
    L2_Z = sumabs2(Z)
    #rank_Z = rank(Z)
    chisq = sumabs2(Y-X*Z)
    println(i, "  ", gamma_list[i], ": ", L0_Z, "  ", L1_Z, " ", L2_Z, " chi^2= ",chisq)# " rank= ", rank_Z, " Z=",Z)
  end
end

function compute_alr_path_stats{T<:Number}(Z_list::Array{T,3}; target=nothing)
  num_it = size(Z_list,1)
  L0_Z_list = zeros(T,num_it)
  L1_Z_list = zeros(T,num_it)
  L2_Z_list = zeros(T,num_it)
  rmsdelta_list = zeros(T,num_it)
  for i in 1:num_it
    Z =  reshape(Z_list[i,:,:],(size(Z_list,2),size(Z_list,3)))
    num_zeros = length(Z)-length(find(Z))
    L0_Z_list[i] = num_nonzeros = length(find(Z))
    L1_Z_list[i] = sumabs(Z)
    L2_Z_list[i] = sumabs2(Z)
    if target != nothing
      rmsdelta_list[i] = sqrt(sumabs2(Z.-target)/length(target))
    end
  end
  return(L0_Z_list,L1_Z_list,L2_Z_list,rmsdelta_list)
end

function compute_alr_path_chisq!{T<:Number}(Z_list::Array{T,3},X::Array{T,2},Y::Array{T,2})
  num_it = size(Z_list,1)
  chisq_list = zeros(T,num_it)
  chisq_refit_list = zeros(T,num_it)
  for i in 1:num_it
    Z =  reshape(Z_list[i,:,:],(size(Z_list,2),size(Z_list,3)))
    chisq_list[i] = sumabs2(Y-X*Z)
    (Zp, chisq_refit) = refit_arp_slr_2d(X,Y,Z)
    Z[i,:,:] = Zp
    chisq_refit_list[i] = chisq_refit
  end
  return (chisq_list, chisq_refit_list)
end

function compute_alr_path_chisq!{T<:Number}(Z_list::Array{T,3}, X::Array{T,2}, Y::Array{T,2}, Xcv::Array{T,2}, Ycv::Array{T,2} )
  num_it = size(Z_list,1)
  chisq_list = zeros(T,num_it)
  chisq_refit_list = zeros(T,num_it)
  chisq_cv_list = zeros(T,num_it)
  for i in 1:num_it
    Z =  reshape(Z_list[i,:,:],(size(Z_list,2),size(Z_list,3)))
    (Zp, chisq_refit, chisq_cv) = refit_arp_slr_2d(X,Y,Z,Xcv,Ycv)
    Z[i,:,:] = Zp
    chisq_refit_list[i] = chisq_refit
    chisq_cv_list[i] = chisq_cv
  end
  return (chisq_list, chisq_refit_list, chisq_cv_list)
end

