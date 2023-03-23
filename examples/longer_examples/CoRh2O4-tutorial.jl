# > ![](https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/sunny_logo.jpg)
# _This is a [tutorial](https://github.com/SunnySuite/SunnyTutorials/tree/main/tutorials) 
#  for the [Sunny](https://github.com/SunnySuite/Sunny.jl/) package, 
#  which enables dynamical simulations of ordered and thermally disordered spins with dipole 
#  and higher order moments._
#
# ## Welcome to a Sunny Tutorial on the Diamond Lattice System CoRh<sub>2</sub>O<sub>4</sub>
# **Script**: Diamond Lattice Finite Temperature Calculation <br>
# **Inspired by**: CoRh<sub>2</sub>O<sub>4</sub> Powder 
# (Ge _et al._ https://doi.org/10.1103/PhysRevB.96.064413) <br>
# **Authors**: Martin Mourigal, David Dalbhom <br>
# **Date**: March 21, 2023  (Sunny 0.4.2) <br>
# **Goal**: This script is to calculate the temperature dependence of the magnon excitations in the 
# spin-3/2 Heisenberg Diamond Antiferromagnet and compare to powder-averaged results obtained for 
# the compound CoRh<sub>2</sub>O<sub>4</sub> <br>

# ---
# ### Loading Packages and Custom Functions 
cd("/Users/mourigal/Dropbox (GaTech)/Group/Shared-Scripts/MourigalLab-Git/Sunny-External/")  #src

# #### Loading Packages 
using Sunny, GLMakie, ProgressMeter, Statistics, Random, Brillouin
import Plots ## Plotting tools as import to avoid conflict with GLMakie
Sunny.offline_viewers() 

# #### Defining Custom functions

# Function:quench! Randomizes spin system, go to a target temperature, let the system relax with 
# nrelax steps of the integrator.
function quench!(sys, integrator; kTtarget, nrelax)
    randomize_spins!(sys);
    integrator.kT = kTtarget;
    prog          = Progress(nrelax; dt=10.0, desc="Quenched and now relaxing: ", color=:green);
    for i in 1:nrelax
        step!(sys, integrator)
        next!(prog)
    end 
end

# Function:dwell! Takes the spin system, go to a target temperature, and dwell with ndwell steps 
# of the Integrator.
function dwell!(sys, integrator; kTtarget, ndwell)
    integrator.kT = kTtarget;
    prog          = Progress(ndwell; dt=10.0, desc="Dwelling: ", color=:green);
    for i in 1:ndwell
        step!(sys, integrator)
        next!(prog)
    end 
end

# Function:anneal! Takes a temperature schedule and cool the system through it, with ndwell steps of 
# the Integrator at each temperature on the schedule. Returns the energy at the end of the dweel for 
# each scheduled temperature.
function anneal!(sys,  integrator;  kTschedule, ndwell)
    nspins = prod(size(sys.dipoles));
    ensys  = zeros(length(kTschedule))        
    prog   = Progress(ndwell*length(kTschedule); dt=10.0, desc="Annealing: ", color=:red);
    for (i, kT) in enumerate(kTschedule)
        for _ in 1:ndwell
            step!(sys, integrator)
            next!(prog)  
        end
        ensys[i] = energy(sys)  
    end
    return ensys/nspins   
end

# Function:sample_sf! Samples a Structure Factor: Instant or Dynamical Structure Factor. 
# The integrator is run ndecorr times before each one of the nsamples is measured.
function sample_sf!(sf, sys, integrator; nsamples, ndecorr)
    prog  = Progress(nsamples*ndecorr; dt=10.0, desc="Sampling SF: ", color=:red);
    for _ in 1:nsamples
        for _ in 1:ndecorr 
            step!(sys, integrator)
            next!(prog)
        end
        add_sample!(sf, sys)    # Accumulate the newly sampled structure factor into `sf`
    end
end
    
# Function:powder_average; Powder Averages a Structure Factor. Works for instant or dynamical 
# Structure Factor
function powder_average(dsf, rs, density; η=0.1, mode=:perp, kwargs...)
    prog   = Progress(length(rs); dt=10., desc="Powder Averaging: ", color=:blue);
    output = zeros(Float64, length(rs), length(ωs(dsf)))
    for (i, r) in enumerate(rs)
        qs = spherical_shell(dsf, r, density)
        if length(qs) == 0
            qs = [[0., 0., 0.]] ## If radius is too small, just look at 0 vector
        end
        vals = intensities(dsf, qs, mode; kwargs...)
        vals = broaden_energy(dsf, vals, (ω,ω₀)->lorentzian(ω-ω₀, η))
        output[i,:] .= mean(vals, dims=1)[1,:]
        next!(prog)
    end
    return output
end

# ---
# ### System Definition for CoRh<sub>2</sub>O<sub>4</sub>

# Define the crystal structure of CoRh$_2$O$_4$  in the conventional cell
xtal    = Crystal("CoRh2O4_#109301.cif";symprec=1e-4)
magxtal = subcrystal(xtal,"Co1")
view_crystal(magxtal,6.0)
print_symmetry_table(magxtal, 4.0)

# Define spin and g-factors
valg = 2.0
valS = 3/2;

# Assign Local Hilbert Space
lhs  = [SpinInfo(1,S=valS;g=valg)]
ffs  = [FormFactor(1,"Co2")];

# Create Spin System and Randomize it
sunmode = :dipole
latsize = (10,10,10)
sys     = System(magxtal, latsize, lhs, sunmode; seed=1)
randomize_spins!(sys)
plot_spins(sys, arrowlength=1.0, linewidth=0.5, arrowsize=1.0)

# Define Exchange Interactions 
scaleJ = 0.63
valJ1  = 1.00*scaleJ
set_exchange!(sys, valJ1, Bond(1, 3, [0, 0, 0]) );

# ---
# ### System thermalization to an ordered, yet finite temperature, state

# Define Langevin Integrator and Initialize it 
Δt0           = 0.05/abs(scaleJ*valS); ## Time steps in Langevin
λ0            = 0.1; ## Langevin damping, usually 0.05 or 0.1 is good.
kT0           = 0.01*abs(scaleJ*valS); ## Initialize at some temperature
integrator    = Langevin(Δt0; λ=λ0, kT=kT0); 

# Thermalization 
# Option 1: Quench the system from infinite temperature to a target temperature. 
# Note: this may lead to a poorly thermalized sample
quench!(sys,integrator;kTtarget=kT0,nrelax=10000);

# Option 2: Anneal (according to a temperature schedule) than dwell once reach base
# Note: starting from very high temperature here 

## kTs = [abs(scaleJ)*valS*100 * 0.9^k for k in 0:100];
## anneal!(sys,integrator;kTschedule=kTs,ndwell=500);
## dwell!(sys,integrator;kTtarget=kTs[end],ndwell=2000);

# Plot the resulting spin system to check ordering in real space
plot_spins(sys, arrowlength=1.0, linewidth=0.5, arrowsize=1.0)

# --- 
# ### Calculation of Neutron Scattering Responses


# #### Fourier transformed instantaneous two-point correlation functions

# Calculate the Instantaneous/Equal-time Structure Factor
@time eqsf = InstantStructureFactor(sys);

# If desired, add additional samples by decorrelating and then re-calculating the eqsf
nsamples   = 0; 
ndecorr    = 1000; 
@time sample_sf!(eqsf, sys, integrator; nsamples=nsamples, ndecorr=ndecorr);

# Project onto a constant Q-Slice in momentum space 
nQpts  = 200;
Qxpts  = range(-10.0, 10.0, length=nQpts);
Qypts  = range(-10.0, 10.0, length=nQpts);
qz     = 1.0;
Qpts   = [[qx, qy, qz] for qx in Qxpts, qy in Qypts];
@time  iq = instant_intensities(eqsf, Qpts, :perp; interpolation = :none, kT=integrator.kT, formfactors=ffs); 

# Plot the resulting I(Q)
fig1=Plots.heatmap(Qxpts, Qypts, iq[:,:]';
    color=cgrad(:viridis, scale = :lin, rev = false),
    clims=(0,maximum(iq)/20), 
    xlabel="Momentum Transfer Qx (r.l.u)", xlabelfontsize=10, 
    ylabel="Momentum Transfer Qy (r.l.u)", ylabelfontsize=10, 
    size=(800,550), fmt=:png,
    margin=5Plots.PlotMeasures.mm, format=:png
);
display(fig1)


# #### Dynamical and energy-integrated two-point correlation functions

# Calculate the Time Traces and Fourier Transform: Dynamical Structure Factor (first sample)
ωmax       = 6.0;  # Maximum  energy to resolve
nω         = 100;  # Number of energies to resolve
@time dsf  = DynamicStructureFactor(sys; Δt=Δt0, nω=nω, ωmax=ωmax, process_trajectory=:symmetrize); 

# If desired, add additional samples by decorrelating and then re-calculating the dsf
nsamples      = 0; # Aditional Samples
ndecorr       = 1000;
@time sample_sf!(dsf, sys, integrator; nsamples=nsamples, ndecorr=ndecorr);

# Can use the Brillouin package for help on determining high symmetry points 
kp        = irrfbz_path(227,[[1,0,0], [0,1,0], [0,0,1]]);
kpc       = cartesianize(kp)

# Project onto a constant QE-Slice in momentum-energy space. 
densQpts  = 50;
symQpts   = [[0.75, 0.75, 0.00],  # List of wave vectors that define a path
            [0.00, 0.00, 0.00],
            [0.50, 0.50, 0.50],
            [0.50, 1.00, 0.00],
            [0.00, 1.00, 0.00],
            [0.25, 1.00, 0.25],
            [0.00, 1.00, 0.00],
            [0.00,-4.00, 0.00]];
(Qpts, symQind) = connected_path(symQpts,densQpts);
@time  iqw = intensities(dsf, Qpts, :perp; interpolation = :none, kT=integrator.kT, formfactors=ffs); 

# If desired, broaden the dsf in energy
η     = 0.02 ## Lorentzian energy broadening parameter
iqwc  = broaden_energy(dsf, iqw, (ω, ω₀) -> lorentzian(ω-ω₀, η));

# If desired, calculated the energy-integrated structure factor
@time  iqt = instant_intensities(dsf, Qpts, :perp; interpolation = :none, kT=integrator.kT, formfactors=ffs); 

# Plot the resulting I(Q,W)    
fig2=Plots.heatmap(1:size(iqwc, 1), ωs(dsf), iqwc';
    color=cgrad(:viridis, scale = :lin, rev = false),
    clims=(0,maximum(iqwc)/20000), 
    xticks = (symQind, symQpts), xrotation=35, xtickfontsize=8,
    xlabel="Momentum Transfer (r.l.u)",
    ylabel="Energy Transfer (meV)", 
    size=(800,550), fmt=:png,
    margin=5Plots.PlotMeasures.mm, format=:png
);
display(fig2)

# Projection into a powder-averaged neutron scattering intensity 
Qmax       = 3.5;
nQpts      = 100;
Qpow       = range(0, Qmax, length=nQpts);  
η          = 0.1;
sphdensity = 4.0;
@time pqw = powder_average(dsf, Qpow, sphdensity; η, kT=integrator.kT, formfactors=ffs);

# Plot resulting Ipow(Q,W)    
fig3=Plots.heatmap(Qpow, ωs(dsf), pqw';
    color=cgrad(:viridis, scale = :log, rev = false),
    clims=(0,maximum(pqw)/100), 
    xlabel="|Q| (Å⁻¹)",
    ylabel="Energy Transfer (meV)", 
    size=(800,550), fmt=:png,
    margin=5Plots.PlotMeasures.mm, format=:png
); 
display(fig3)

# --- 
# ### Calculation of temperature-dependent powder average spectrum

# Define a temperature schedule
kTs        = [60 40 25 20 15 12 10 4]/11.602;
pqw_res    = Array{Matrix{Float64}}(undef, length(kTs)); 
iqw_res    = Array{Matrix{Float64}}(undef, length(kTs)); 
for (i, kT) in enumerate(kTs)
    dwell!(sys,integrator;kTtarget=kT,ndwell=1000);
    dsf_loc     = DynamicStructureFactor(sys; Δt=2Δt0, nω=nω, ωmax=ωmax, process_trajectory=:symmetrize); 
    iqw_res[i]  = intensities(dsf_loc, Qpts, :perp; interpolation = :none, kT=kT, formfactors=ffs); 
    pqw_res[i]  = powder_average(dsf_loc, Qpow, sphdensity; η, kT=kT, formfactors=ffs);
end

# Plot the resulting Ipow(Q,W) as a function of temperature,
# to compare with Fig.6 of https://arxiv.org/abs/1706.05881
for i = 8:-1:1
    fig4 = Plots.heatmap(Qpow, ωs(dsf), pqw_res[i]';
        color=cgrad(:viridis, scale = :lin, rev = false),
        clims=(0,10), 
        xlabel="|Q| (Å⁻¹)",
        ylabel="Energy Transfer (meV)", 
        size=(800,550), fmt=:png,
        margin=5Plots.PlotMeasures.mm, format=:png
    );
    display(fig4)
end

