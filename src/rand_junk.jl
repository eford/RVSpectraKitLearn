
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

doppler_comp = calc_doppler_component(lambda,obs, sigmasq_obs=1/150000.0*ones(length(lambda)), rho=0.1)
genpca_out = fit_gen_pca_eford(obs,doppler_comp)

doppler_quadratic_term = calc_doppler_quadratic_term(lambda,obs, sigmasq_obs=1/150000.0*ones(length(lambda)), rho=0.1)
gen2pca_out = fit_gen2_pca_eford(obs,doppler_comp,doppler_quadratic_term)








