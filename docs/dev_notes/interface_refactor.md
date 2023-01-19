
## Sampling

Need to figure out a new interface for samplers and Monte Carlo.

## Inhomogeneities

After discussion with David, we are planning this interface:

```julia
set_exchange!(sys, J, i) # as usual
get_exchange(sys, i) # for convenience, allow to inspect 3x3 matrices

set_external_field!(sys, h0)
set_external_field_at!(sys, h, idx) # this is easy to support, keeps h0 as the background
get_external_field_at!(sys, idx)    # might as well have this too

# New function, enter "inhomogeneous mode". Can't go back.
enable_inhomogeneous_exchange!(sys) 

set_exchange_at!(sys, J, idx)  # this becomes possible
get_exchange_at(sys, idx)      # this too

set_exchange!(sys, J, i) # ERROR! can no longer define symmetry-propagated exchange
get_exchange(sys, i)     # ERROR! I think in inhomogeneous mode this becomes confusing.
get_reference_exchange_at(sys, idx) # Let users see original exchange matrices prior to inhomogeneous modifications?

set_external_field!(sys, h1) # still allowed; overwrites external field at every site
```
