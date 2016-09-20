function get_filenames(; num_max_files::Integer = 0)
  files = readdir()
  datafiles = String[]
  strengths = Float64[]
  snrs = Float64[]
  phases = Float64[]
  for f in files
    m = match(r"(\w+)_(\d+)\w+_(\d+)\w+_(-?\d+\.\d*)\w+_(-?\d+\.\d*)\.csv",f)
    if m==nothing continue end
    effect_name = m.captures[1]
    push!(strengths,parse(m.captures[2]))
	  push!(snrs,parse(m.captures[3]))
    sampling = parse(m.captures[4])
    push!(phases,parse(m.captures[5]))
    push!(datafiles,f)
    if num_max_files>0 && length(datafiles)==num_max_files break end
  end
  @assert(length(datafiles)>0)
  idx = sortperm(phases)
  return (datafiles[idx],phases[idx])
end

function read_filelist(datafiles::Vector{String})
  data = readcsv(datafiles[1],skipstart=1);
  lambda = data[:,1]
  p = size(data,1)
  n = length(datafiles)
  X = zeros(p,n);
  for i in 1:n
    X[:,i] = readcsv(datafiles[i],skipstart=1)[:,2];
  end
  return (lambda,X)
end

# Mirrors code in  Sept 8, 2016 email from Allen Davis
function make_noisy_spectrum{T<:Real}(spec::Array{T,1}, snr::T; sampling::T = 1.0)
  snr_per_pixel = float(snr) / sqrt(sampling)
  orig_mean = mean(spec)
  scaled_spec = spec*snr*snr/orig_mean
  noisy_spec = scaled_spec + randn(length(spec)).*sqrt(scaled_spec) #  np.random.normal(scaled_spec,np.sqrt(scaled_spec))
  return orig_mean*noisy_spec/(snr*snr)
end

function make_noisy_spectra{T<:Real}(spec::Array{T,2}, snr::T; sampling::T = 1.0)
  output = similar(spec)
  for i in 1:size(spec,2)
    output[:,i] = make_noisy_spectrum(spec[:,i], snr, sampling=sampling)
  end
  return output
end






