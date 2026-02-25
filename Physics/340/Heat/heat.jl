module Heat
using Symbolics
using DomainSets


# Ch. 18

#Expansion
@symbols α β ΔL ΔV ΔT L0 V0
linear_expansion = ΔL ~ L0*α*ΔT
volumetric_expansion = ΔV ~ V0*β*ΔT

# Heat Capacity & Specific Heat
@symbols Q C Tf Ti c m
heat_capacity = Q ~ C * (Tf - Ti)
specific_heat = Q ~ c * m * (Tf - Ti)

# Heat of Transformation
@symbols Q L m
heat_of_trans = Q ~ L * m

#Work, ΔVol.
@symbols W V Vi Vf p
work_volumetric = W ~ Integral(V in ClosedInterval(Vi,Vf))(p)

# First Law of Thermo
@symbols ΔU(..) δQ(..) δW(..)
first_law_thermo = ΔU ~ δQ - δW
#=
Adiabatic Processes:        Q = 0, ΔU = -W
Constant-volume Processes:  W = 0, ΔU = Q
Cyclical Processes:         ΔU = 0, Q = W
Free Expansion:             Q = W = ΔU = 0
=#

#Conduction, Convection, & Radiation
@symbols Pcond Prad Pabs Q t k A TH TC L σ ε T Tenv
conduction_short    = Pcond ~  Q/t
conduction_full     = Pcond ~ k*A*((TH-TC)/L)

radiation_heat      = σ*ε*A*T^4
absorbtion_heat     = σ*ε*A*Tenv^4

# Chap 19

# Ideal Gases
@symbols p V n R T N k W Vf Vi ΔV
ideal_gas_law_moles     = p*V ~ n*R*T
ideal_gas_law_molecules = p*V ~ N*k*T
ideal_gas_work          = W ~ n*R*T*log(Vf/Vi)
