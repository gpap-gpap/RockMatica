(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["RockMatica`MultiFluid`"]


(* ::Input::Initialization:: *)
rpFluidMix::usage = "rpFluidMix[K1, K2, s1, capillaryParameter->q]"; 
rpViscosityMix::usage = "rpViscosityMix[\[Eta]1, \[Eta]2, s1, capillaryParameter->q, \
BrooksCoreyParameters->{\[Lambda],{sw_red, snw_red}}]"; 
rpBrooksCorey::usage = "rpBrooksCorey[\[Lambda], {sw_red, snw_red}]"; 
rpTimeConstantMix::usage = "rpTimeConstantMix[\[Eta]1, \[Eta]2, s1, capillaryParameter->q, \
BrooksCoreyParameters->{\[Lambda],{sw_red, snw_red}}]"; 
rpDensityMix::usage = "rpDensityMix[\[Rho]1, \[Rho]2, s1]"


(* ::Input::Initialization:: *)
Begin["Private`"]


numerical = Chop@N@#&;


(* ::Input::Initialization:: *)
Options[rpFluidModuliMix] = {capillaryParameter -> 1}; 
rpFluidModuliMix[K1_, K2_, s1_,opts:OptionsPattern[]] := 
   Module[{q = OptionValue[capillaryParameter]}, 
  (s1 + q*(1 - s1))/(s1/K1 + q*((1 - s1)/K2))//numerical]; 
rpBrooksCorey[\[Lambda]_, {swr_,snwr_}] := Module[{kw, knw, seff}, 
    seff = Clip[(#1 - swr)/(1 - swr - snwr),{0, 1}] & ; 
     kw = Clip[seff[#1]^((2 + 3*\[Lambda])/\[Lambda]),{0, 1}]&; 
  knw = Clip[(1 - seff[#1])^2*(1 - seff[#1]^((2 + \[Lambda])/\[Lambda])), {0, 1}]&; 
     {kw[#1], knw[#1]}&]; 
Options[rpViscosityMix] ={BrooksCoreyParameters ->{1, {0, 0}}}; 
rpViscosityMix[\[Eta]1_, \[Eta]2_, s1_,opts:OptionsPattern[{rpViscosityMix,rpFluidModuliMix}]] := 
   Module[{q = OptionValue[capillaryParameter], 
     k = rpBrooksCorey@@OptionValue[BrooksCoreyParameters]}, (s1 + q*(1 - s1))/(k[s1][[1]]/\[Eta]1 + q*(k[s1][[2]]/\[Eta]2))//numerical]; 
Options[rpTimeConstantMix] = 
   {capillaryParameter -> 1,BrooksCoreyParameters -> {1, {0, 0}}}; 
rpTimeConstantMix[\[Eta]1_, \[Eta]2_,s1_, opts:OptionsPattern[
      {rpViscosityMix,rpTimeConstantMix}]] := 
   Module[{q = OptionValue[capillaryParameter], 
     bk = OptionValue[BrooksCoreyParameters]}, 
    (1/\[Eta]1)*rpViscosityMix[\[Eta]1,\[Eta]2, s1, capillaryParameter -> q,BrooksCoreyParameters ->bk//numerical]]; 
rpDensityMix[\[Rho]1_, \[Rho]2_, s1_]:=\[Rho]1 s1 + (1-s1)\[Rho]2;


Options[rpMultiFluid] = {FluidSaturation->0., PatchParameter->1}
rpMultiFluid[a_RockMatica`Base`rpFluid, b_RockMatica`Base`rpFluid, OptionsPattern[]]:=
	Module[{bulk, visc, dens},
	
	rpMultiFluid[fluidName]:=<|
	"FluidName"->fluidName,
	"FluidModulus"->fluidMod, 
	"FluidViscosity"->fluidVisc,
	"FluidDensity"->fluidDens|>
	]


(* ::Input::Initialization:: *)
End[]


(* ::Input::Initialization:: *)
EndPackage[]
