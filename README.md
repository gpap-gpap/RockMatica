# RockMatica
A rock physics package for Mathematica. As it stands two isotropic and TI anisotropic elastic tensors are supported. Example workflow:

## Define rocks/fluids
The command 

```Mathematica 
<<RockMatica` (* Assumes the repository is cloned in $UserBaseDirectory/RockMatica *)
rpRock["RockName", rock properties]
```
defines an `npRock` object (essentially an association with ``` "ParameterName"->ParameterValue ```). The rock properties supported can be show by running ``` ?rpRock```. Essentially, if ran with different number of parameters the program assumes that the rock defined is isotropic (if only bulk and shear modulus are given) or TI with symmetry along the z-axis (HTI).

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

The option "MicrocrackParameters" refers to the viscoelastic rock physics model based on Eshelby inclusions detailed in Chapman, Mark, Sergei V. Zatsepin, and Stuart Crampin. "Derivation of a microstructural poroelastic model." Geophysical Journal International 151.2 (2002): 427-451.

This way, a new rock physics model can be implemented at the rock definition stage and calculation of its output can be performed by checking if the relevant options are defined.


## Calculate Elastic Tensor


# To-Do
- count length of rpRock input to ensure that the rock is iso/aniso tropic. OR add a flag ```Mathematica anisotropic->True```
- add multi fluid description
