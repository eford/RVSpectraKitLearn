function get_filenames()
  files = readdir()
  datafiles = String[]
  strengths = Float64[]
  snrs = Float64[]
  phases = Float64[]
  for f in files
    m = match(r"(\w+)_(\d+)\w+_(\d+)\w+_(-?\d+\.\d*)\.csv",f)
    if m==nothing continue end
    effect_name = m.captures[1]
    push!(strengths,parse(m.captures[2]))
	push!(snrs,parse(m.captures[3]))
    push!(phases,parse(m.captures[4]))
    push!(datafiles,f)
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






