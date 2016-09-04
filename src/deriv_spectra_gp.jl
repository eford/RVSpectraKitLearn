# WARNING: Work in progress
# E.g., probably needs to update function names/signatures
# Need way to optimize covar to be useful

using PDMats

function kernel_matern32(d::Float64; rho::Float64 = 1.0, sigmasq::Float64 = 1.0)
  x = abs(sqrt(3)*d/rho)
  sigmasq * (1+x) * exp(-x)
end

function dkerneldx_matern32(d::Float64; rho::Float64 = 1.0, sigmasq::Float64 = 1.0)
  x = abs(sqrt(3)*d/rho)
  -sigmasq * (-3*d/rho^2) * exp(-x)
end

function d2kerneldx2_matern32(d::Float64; rho::Float64 = 1.0, sigmasq::Float64 = 1.0)
  x = abs(sqrt(3)*d/rho)
  -sigmasq * (3/rho^2) * (1-x)*exp(-x)
end

function make_kernel_data{T}(x::AbstractArray{T,1}, kernel::Function;
			sigmasq_obs::AbstractArray{T,1} = zeros(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)
  @assert length(x) == length(sigmasq_obs)
  K = diagm(sigmasq_obs)
  for i in 1:size(K,1)
      for j in 1:size(K,2)
	      K[i,j] += kernel(x[i]-x[j], sigmasq=sigmasq_cor, rho=rho)
	  end
  end
  return PDMat(K)
end

function make_kernel_obs_pred{T}(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1}, kernel::Function;
			sigmasq_cor::T = 1.0, rho::T = 1.0)
  K = zeros(length(xobs),length(xpred))
  for i in 1:length(xobs)
      for j in 1:length(xpred)
	      K[i,j] += kernel(xobs[i]-xpred[j], sigmasq=sigmasq_cor, rho=rho)
	  end
  end
  return K
end

function make_kernel_matern32_data{T}(x::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = zeros(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)
  make_kernel_data(x, kernel_matern32, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_matern32_obs_pred{T}(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)
  make_kernel_obs_pred(xobs, xpred, kernel_matern32, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_dmatern32dx_obs_pred{T}(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)
  make_kernel_obs_pred(xobs, xpred, dkerneldx_matern32, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_d2matern32dx2_obs_pred{T}(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)
  make_kernel_obs_pred(xobs, xpred, d2kerneldx2_matern32, sigmasq_cor=sigmasq_cor, rho=rho)
end

function predict_mean{T}(xobs::AbstractArray{T,1}, yobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(xobs)),	sigmasq_cor::T = 1.0, rho::T = 1.0)
  kobs = make_kernel_matern32_data(xobs, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  kobs_pred = make_kernel_matern32_obs_pred(xobs,xpred, sigmasq_cor=sigmasq_cor, rho=rho)
  alpha = kobs \ yobs
  pred_mean = kobs_pred' * alpha
end

function predict_deriv{T}(xobs::AbstractArray{T,1}, yobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(xobs)),	sigmasq_cor::T = 1.0, rho::T = 1.0)
  kobs = make_kernel_matern32_data(xobs, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  kobs_pred_deriv = make_kernel_dmatern32dx_obs_pred(xobs,xpred, sigmasq_cor=sigmasq_cor, rho=rho)
  alpha = kobs \ yobs
  pred_deriv = kobs_pred_deriv' * alpha
end

function predict_deriv2{T}(xobs::AbstractArray{T,1}, yobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(xobs)),	sigmasq_cor::T = 1.0, rho::T = 1.0)
  kobs = make_kernel_matern32_data(xobs, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  kobs_pred_deriv2 = make_kernel_d2matern32dx2_obs_pred(xobs,xpred, sigmasq_cor=sigmasq_cor, rho=rho)
  alpha = kobs \ yobs
  pred_deriv = kobs_pred_deriv2' * alpha
end

function predict_mean_and_derivs{T}(xobs::AbstractArray{T,1}, yobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(xobs)),	sigmasq_cor::T = 1.0, rho::T = 1.0)
  kobs = make_kernel_matern32_data(xobs, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  alpha = kobs \ yobs
  kobs_pred = make_kernel_matern32_obs_pred(xobs,xpred, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_mean = kobs_pred' * alpha
  kobs_pred_deriv = make_kernel_dmatern32dx_obs_pred(xobs,xpred, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_deriv = kobs_pred_deriv' * alpha
  kobs_pred_deriv2 = make_kernel_d2matern32dx2_obs_pred(xobs,xpred, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_deriv2 = kobs_pred_deriv2' * alpha
  return (pred_mean, pred_deriv, pred_deriv2)
end

function gp_marginal{T}(xobs::AbstractArray{T,1}, yobs::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(xobs)),	sigmasq_cor::T = 1.0, rho::T = 1.0)
  kobs = make_kernel_matern32_data(xobs, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  -0.5*( invquad(kobs, yobs) + logdet(kobs) + length(xobs)*log(2pi) )
end

function calc_gp_marginal_on_segments{T}(lambda::AbstractArray{T,1}, flux::AbstractArray{T,1};
                                sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(lambda)),	sigmasq_cor::T = 1.0, rho::T = 1.0,
                                half_chunck_size::Integer = 100)
  @assert length(lambda) == length(flux) == length(sigmasq_obs)
  println("# sigmasq_obs[1,1] = ", sigmasq_obs[1,1], " sigmasq_cor= ", sigmasq_cor, " rho= ", rho)
  output = 0.0
  num_seg = convert(Int64,ceil(length(lambda)/half_chunck_size)-1)
  for i in 1:num_seg
    idx_begin = 1+half_chunck_size*(i-1)
    idx_end = min(half_chunck_size*(i+1), length(lambda))
    write_idx_begin = idx_begin + div(half_chunck_size,2)
    write_idx_end = idx_end - div(half_chunck_size,2)
    if i==1 write_idx_begin=1 end
    if i==num_seg
      idx_begin = max(1,idx_end-2*half_chunck_size)
      write_idx_end=length(lambda)
    end
    #println("# i= ",i,": ", idx_begin, " - ", idx_end, " -> ", write_idx_begin, " - ", write_idx_end)
    #output[write_idx_begin:write_idx_end] = predict_gp(view(lambda,idx_begin:idx_end), view(flux,idx_begin:idx_end), view(lambda,write_idx_begin:write_idx_end), sigmasq_obs=view(sigmasq_obs,idx_begin:idx_end), sigmasq_cor=sigmasq_cor, rho=rho)
    output += gp_marginal(lambda[idx_begin:idx_end], flux[idx_begin:idx_end], sigmasq_obs=sigmasq_obs[idx_begin:idx_end], sigmasq_cor=sigmasq_cor, rho=rho)
  end
  output
end

function gp_marginal_wrapper(param::Vector)
   @assert length(param) == 2
   idx_min = 40000
   idx_max = 45000
   println("# wrapper: ", param)
   calc_gp_marginal_on_segments(lambda[idx_min:idx_max],vec(mean(obs[idx_min:idx_max],2)), sigmasq_obs = ones(length(lambda[idx_min:idx_max]))/150000,
                    sigmasq_cor = exp(param[1]), rho = exp(param[2]) )
end





function calc_doppler_component{T}(lambda::AbstractArray{T,1}, flux::AbstractArray{T,1};
                                sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(lambda)),	sigmasq_cor::T = 1.0, rho::T = 1.0,
                                half_chunck_size::Integer = 100)
   lambda.*calc_gp_on_segments(predict_deriv,lambda, flux, sigmasq_obs=sigmasq_obs,	sigmasq_cor=sigmasq_cor, rho=rho, half_chunck_size=half_chunck_size)
end

function calc_gp_on_segments{T}(predict_gp::Function, lambda::AbstractArray{T,1}, flux::AbstractArray{T,1};
                                sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(lambda)),	sigmasq_cor::T = 1.0, rho::T = 1.0,
                                half_chunck_size::Integer = 100)
  @assert length(lambda) == length(flux) == length(sigmasq_obs)
  output = Array(Float64,length(lambda))
  num_seg = convert(Int64,ceil(length(lambda)/half_chunck_size)-1)
  for i in 1:num_seg
    idx_begin = 1+half_chunck_size*(i-1)
    idx_end = min(half_chunck_size*(i+1), length(lambda))
    write_idx_begin = idx_begin + div(half_chunck_size,2)
    write_idx_end = idx_end - div(half_chunck_size,2)
    if i==1 write_idx_begin=1 end
    if i==num_seg
      idx_begin = max(1,idx_end-2*half_chunck_size)
      write_idx_end=length(lambda)
    end
    #println("# i= ",i,": ", idx_begin, " - ", idx_end, " -> ", write_idx_begin, " - ", write_idx_end)
    #output[write_idx_begin:write_idx_end] = predict_gp(view(lambda,idx_begin:idx_end), view(flux,idx_begin:idx_end), view(lambda,write_idx_begin:write_idx_end), sigmasq_obs=view(sigmasq_obs,idx_begin:idx_end), sigmasq_cor=sigmasq_cor, rho=rho)
    output[write_idx_begin:write_idx_end] = predict_gp(lambda[idx_begin:idx_end], flux[idx_begin:idx_end], lambda[write_idx_begin:write_idx_end], sigmasq_obs=sigmasq_obs[idx_begin:idx_end], sigmasq_cor=sigmasq_cor, rho=rho)
  end
  output
end

function calc_doppler_component{T}(lambda::AbstractArray{T,1}, flux::AbstractArray{T,2};
                                sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(xobs)),	sigmasq_cor::T = 1.0, rho::T = 1.0)
  doppler_basis = calc_doppler_component(lambda,vec(mean(flux,2)),sigmasq_obs=sigmasq_obs,sigmasq_cor=sigmasq_cor,rho=rho)
end

function calc_doppler_quadratic_term{T}(lambda::AbstractArray{T,1}, flux::AbstractArray{T,1};
                                sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(lambda)),	sigmasq_cor::T = 1.0, rho::T = 1.0,
                                half_chunck_size::Integer = 100)
   0.5*lambda.^2.*calc_gp_on_segments(predict_deriv2,lambda, flux, sigmasq_obs=sigmasq_obs,	sigmasq_cor=sigmasq_cor, rho=rho, half_chunck_size=half_chunck_size)
end

function calc_doppler_quadratic_term{T}(lambda::AbstractArray{T,1}, flux::AbstractArray{T,2};
                                     sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(xobs)),	sigmasq_cor::T = 1.0, rho::T = 1.0)
  doppler_quad_term = calc_doppler_quadratic_term(lambda,vec(mean(flux,2)),sigmasq_obs=sigmasq_obs,sigmasq_cor=sigmasq_cor,rho=rho)
end


function test_gp_deriv()
  #lambda_small = collect(linspace(1000.0,1100.0,101))
  #yobs_small = sin(2pi*lambda/40)+ 0.3*randn(length(lambda_small))
  #lambda_pred = collect(linspace(1000.0,1100.0,201))
  lambda_small = lambda[10000:10100]
  yobs_small = obs[10000:10100,1]
  lambda_pred = lambda_small
  (f,df,df2) = predict_mean_and_derivs(lambda_small,yobs_small,lambda_pred, sigmasq_obs=1.0/150000.0*ones(length(lambda_small)), sigmasq_cor=0.1,rho = 0.01)
  plot(lambda_small,(yobs_small),"r.")
  plot(lambda_pred,f,"b-")
  plot(lambda_pred,df,"g-")
  #plot(lambda_pred,df2,"m-")
end

#=
test_gp_deriv()

doppler_comp = calc_doppler_component(lambda,obs, sigmasq_obs=1/150000.0*ones(length(lambda)), rho=0.1)
genpca_out = fit_gen_pca_eford(obs,doppler_comp)

doppler_quadratic_term = calc_doppler_quadratic_term(lambda,obs, sigmasq_obs=1/150000.0*ones(length(lambda)), rho=0.1)
genpca_quad_out = fit_gen_pca_eford(obs,doppler_comp,qaudratic_term=doppler_quadratic_term)
=#


