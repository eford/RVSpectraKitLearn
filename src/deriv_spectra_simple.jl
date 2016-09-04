# Functions to estimate the derivative(s) of the mean spectrum (probably could do better w/ spline or GP)
function calc_deriv_simple{T}(x::AbstractArray{T,1})
  @assert length(x)>=3
  dx = similar(x)
  dx[1] = x[2]-x[1]
  dx[end] = x[end]-x[end-1]
  for i in 2:(length(x)-1)
     dx[i] = (x[i+1]-x[i-1])/2
  end
  return dx
end

function calc_deriv2_simple{T}(x::AbstractArray{T,1})
  @assert length(x)>=3
  dx2 = similar(x)
  dx2[1] = x[3]-2*x[2]+x[1]
  dx2[end] = x[end]-2*x[end-1]+x[end-2]
  for i in 2:(length(x)-1)
     dx2[i] = (x[i+1]+x[i-1]-2*x[i])
  end
  return dx2
end

function calc_doppler_component_simple{T}(lambda::AbstractArray{T,1}, flux::AbstractArray{T,1})
  @assert length(lambda) == length(flux)
  dlambdadpix = calc_deriv_simple(lambda);
  dfluxdpix = calc_deriv_simple(flux);
  doppler_basis = dfluxdpix.*(lambda./dlambdadpix)
end

function calc_doppler_component_simple{T}(lambda::AbstractArray{T,1}, flux::AbstractArray{T,2})
  doppler_basis = calc_doppler_component_simple(lambda,vec(mean(flux,2)))
end

function calc_doppler_quadratic_term_simple{T}(lambda::AbstractArray{T,1}, flux::AbstractArray{T,1})
  @assert length(lambda) == length(flux)
  dlambdadpix = calc_deriv_simple(lambda);
  d2fluxdpix2 = calc_deriv2_simple(flux);
  doppler_quad_term = 0.5*d2fluxdpix2.*(lambda.^2./dlambdadpix.^2)
end

function calc_doppler_quadratic_term_simple{T}(lambda::AbstractArray{T,1}, flux::AbstractArray{T,2})
  doppler_quad_term = calc_doppler_quadratic_term_simple(lambda,vec(mean(flux,2)))
end


