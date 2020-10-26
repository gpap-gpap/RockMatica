(* ::Package:: *)

(*  Mathematica Package  *)
(* :Title:   RockMatica  *)
(* :Author:  Giorgos Papageorgiou <gpap.gpap@gmail.com> *)
(* :Context: RockMatica`Multifluid` *)
(* :Date:    01/11/2019  *)

(* :Version: 12 *)
(* :Copyright: (c) 2019 Giorgos Papageorgiou *)


(* ::Input::Initialization:: *)
BeginPackage["RockMatica`MultiFluid`"]


(* ::Input::Initialization:: *)
rpFluidMix::usage="rpFluidMix[fluid1_rpFluid, fluid2_rpFluid] defines a FluidMix object assuming first fluid is wetting and second fluid nonwetting. The fluid mix is an association with the additional options 'BrooksCoreyParameters', 'FluidSaturation; and 'PatchParameter'. The name of the fluid mix follows the convention fluid1Name-fluid2Name-Mix ";


(* ::Input::Initialization:: *)
Begin["Private`"]


Options[rpFluidMix] = {"FluidSaturation"->0., "PatchParameter"->1, "BrooksCoreyParameters"->{1,{0,0}}}
rpFluidMix[a_RockMatica`Base`rpFluid, b_RockMatica`Base`rpFluid, OptionsPattern[]]:=
	Module[{bulk, visc, fluidName, brooksCorey, relperm, dens, s, q, l, K1, K2, eta1, eta2, rho1, rho2},
		fluidName=a["FluidName"]<>"-"<>b["FluidName"]<>"-Mix";
		s = OptionValue["FluidSaturation"];
		q = OptionValue["PatchParameter"];
		l = OptionValue["BrooksCoreyParameters"];
		K1=a["FluidModulus"];
		K2=b["FluidModulus"];
		eta1=a["FluidViscosity"];
		eta2=b["FluidViscosity"];
		rho1=a["FluidDensity"];
		rho2=b["FluidDensity"];
		brooksCorey[lambda_, {swr_,snwr_}] := Module[{kw, knw, seff},
			seff = Clip[(#1 - swr)/(1 - swr - snwr),{0, 1}] &;
			kw = Clip[seff[#1]^((2 + 3*lambda)/lambda),{0, 1}]&;
			knw = Clip[(1 - seff[#1])^2*(1 - seff[#1]^((2 + lambda)/lambda)), {0, 1}]&;
			{kw[#1], knw[#1]}&];
		relperm = (brooksCorey@@l)[s];
		bulk=(s + q*(1 - s))/(s/K1 + q*((1 - s)/K2));
		visc=(s + q*(1 - s))/(relperm[[1]]/eta1 + q*(relperm[[2]]/eta2));
		dens= rho1 s + (1-s)rho2;
		rpFluidMix[fluidName]=<|
			"FluidName"->fluidName,
			"FluidModulus"->bulk, 
			"FluidViscosity"->visc,
			"FluidDensity"->dens,
			"RelativePermeability"->(brooksCorey@@l)[s],
			"FirstFluidSaturation"->s,
			"CapillaryParameter"->q
			|>//N
	]


(* ::Input::Initialization:: *)
End[]


(* ::Input::Initialization:: *)
EndPackage[]
