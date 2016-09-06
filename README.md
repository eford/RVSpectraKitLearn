# RvSpectraKitLearn.jl
======================

[![Build Status](https://travis-ci.org/eford/RvSpectraKitLearn.jl.svg?branch=master)](https://travis-ci.org/eford/RvSpectraKitLearn.jl)

[![Coverage Status](https://coveralls.io/repos/eford/RvSpectraKitLearn.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/eford/RvSpectraKitLearn.jl?branch=master)

[![codecov.io](http://codecov.io/github/eford/RvSpectraKitLearn.jl/coverage.svg?branch=master)](http://codecov.io/github/eford/RvSpectraKitLearn.jl?branch=master)


I implemented the generalized PCA algorithm that I had previously described only conceputally for our proposal.
The code is at https://github.com/eford/RVSpectraKitLearn.jl

Some instructions below for trying it out...

Installation
------------
Download and install a v0.5 release candidate of Julia from http://julialang.org/downloads/

Start julia

```shell
julia
```

Install the package
```julia
Pkg.clone("git@github.com:eford/RvSpectraKitLearn.jl.git")  # Install package
```

Edit the file path_to_spectra.jl to point to the directory with your spectra, perhaps like
```julia
path_to_spectra = "C:\\Users\\eford\\Box Sync\\SOAP simulations\\Aug2016_workshop\\SOAP_Spectra\\planet_10ms_150k"
path_to_spectra = joinpath(homedir(),"SOAP_Spectra/planet_10ms_150k")
```

Test the package
----------------
```julia
include(joinpath(Pkg.dir("RvSpectraKitLearn"),"test","runtests.jl"))
```

Example usage
-------------
You can see how to use as an example (very similar to tests), either from julia
```julia
include(joinpath(Pkg.dir("RvSpectraKitLearn"),"examples","gpca_ex1.jl"))
```
or from the shell
```shell
julia examples/gpca_ex1.jl
```

Interoperability
-----------------
If you want to call Python from Julia, see https://github.com/stevengj/PyCall.jl .
If you want to call Julia from python, see https://github.com/JuliaInterop/pyjulia  (I haven't tested this).


