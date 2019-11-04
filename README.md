# RockMatica
A rock physics package for Mathematica. As it stands two isotropic and TI anisotropic elastic tensors are supported. Example workflow:

## Define rocks/fluids
The command 

```Mathematica 
<<RockMatica` (* Assumes the repository is cloned in $UserBaseDirectory/RockMatica *)
rpRock["RockName", rock properties]
```
defines a rock object (essentially an association with "ParameterName"->ParameterValue).
```Mathematica 
rpFluid["FluidName", fluid properties]
```

## Calculate Elastic Tensor


# To Do-s
- count length of rpRock input to ensure that the rock is iso/aniso tropic. OR add a flag ```Mathematica anisotropic->True```
