# RockMatica
A rock physics and wave modelling package for Mathematica. As it stands isotropic and TI anisotropic elastic and viscoelastic tensors are supported. The package consists of subpackages:

- `Base`: base function definitions for rocks and fluids,
- `Moduli`: definitions of elastic and viscoelastic moduli
- `MultiFluid`: definitions of properties of effective fluids from mixtures
- `Waves`: definitions of wavefield modelling code using the reflectivity method
- `Wedge`: definitions of convolution modelling code and inversion using pursuit methods

Example workflow:

## Define rocks/fluids
The command 

```Mathematica 
<<RockMatica` (* Assumes the repository is cloned in $UserBaseDirectory/RockMatica *)
rpRock["RockName", rock properties]
```
defines an `rpRock` object (an association with ``` "ParameterName"->ParameterValue ```). The rock properties supported can be show by running ``` ?rpRock```. Depending on the number of parameters given the program assumes that the rock defined is isotropic (if only bulk and shear modulus are given) or TI with symmetry along the z-axis (VTI).

Likewise, to define a fluid, one calls
```Mathematica 
rpFluid["FluidName", fluid properties]
```
The associations defined with these two operations, can be recalled using only the rock(resp. fluid) name. So after running
```Mathematica
rpRock["sandstone", 12.4, 11, 36, .26, 2.6, MicrocrackParameters -> {.013, 10^-5, 1}]
```
one can recall the rock named "sandstone" by simply running

```Mathematica
rpRock["sandstone"]
Out: <|"RockName" -> "sandstone", "DryModulus" -> 12.4, "ShearModulus" -> 11., "MineralModulus" -> 36., "Porosity" -> 0.26, "MineralDensity" -> 2.6, 
 "MicrocrackParameters" -> <|"CrackDensity" -> 0.013,  "AspectRatio" -> 0.00001, "ReferenceFrequency" -> 1.|>, 
 "FractureParameters" -> <|"FractureDensity" -> None, "FractureFrequencyRatio" -> None|>|>
```

The reason rocks and fluids are defined not as associations but with numerical values is to reduce the chance of mis-spelling seeing as the models rely on spelling of association keys to be uniform. The option "MicrocrackParameters" refers to the viscoelastic rock physics model based on Eshelby inclusions detailed in Chapman, Mark, Sergei V. Zatsepin, and Stuart Crampin. "Derivation of a microstructural poroelastic model." Geophysical Journal International 151.2 (2002): 427-451.

This way, a new rock physics model can be implemented at the rock definition stage and calculation of its output can be performed by checking if the relevant options are defined.


## Calculate Elastic Tensor

## Create a Wedge Model and Invert Using Basis Pursuit


# To-Do
- add error messaging
