(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["RockMatica`AuxiliaryFunctions`"]


(* ::Input::Initialization:: *)
ClearAll[rpVelocity, rpAttenuation, rpPlotGrid];
rpVelocity::usage=""; 
rpAttenuation::usage = ""; 
rpPlotGrid::usage = "rpPlotGrid[PlotMatrix, sizeX, sizeY] returns an aligned plot combining the plots in the plot matrix. Works best with framed plots and fails to align Graphics objects with plots" ; 


(* ::Input::Initialization:: *)
Begin["Private`"]


(* ::Input::Initialization:: *)
Options[rpVelocity] = 
   {solidDensity -> 2650, 
    fluidDensity -> 1000}; 
SetAttributes[rpVelocity,HoldFirst];
rpVelocity[a:(_RockPhysics`Moduli`rpSquirtIsotropic|_RockPhysics`Moduli`rpGassmannIsotropic),opts:OptionsPattern[]]:=Module[{\[Phi] = {Sequence@@Unevaluated[a]}[[4]],\[Rho]}, 
\[Rho]=(1 - \[Phi])*OptionValue[solidDensity] + \[Phi]*OptionValue[fluidDensity];
    Re[{(Sqrt[(#1 + (4/3)*#2)/\[Rho]]),Sqrt[#1/\[Rho]]}&@@(Evaluate@a)]
];
rpVelocity[a:{k_,m_},\[Phi]_, opts:OptionsPattern[]] :=Module[{\[Rho]}, 
\[Rho]=(1 - \[Phi])*OptionValue[solidDensity] + \[Phi]*OptionValue[fluidDensity];
    Re[{(Sqrt[(k + (4/3)*m)/\[Rho]]),Sqrt[m/\[Rho]]}]];
rpVelocity[_]:=$Failed;
rpAttenuation[a:{k_,m_}]:= 
   Module[{Pmod=(k+(4/3)*m &)@@a}, 
    {Im[Pmod]/Re[Pmod],Im@a[[2]]/Re@a[[2]]}];
Options[plotGrid] = 
   {ImagePadding -> 50}; 
rpPlotGrid[l_List, w_, h_, 
    opts:OptionsPattern[]] := 
   Module[{nx, ny, sidePadding = 
      OptionValue[plotGrid, 
       ImagePadding], 
     topPadding = 0, widths, 
     heights, dimensions, 
     positions, frameOptions = 
      FilterRules[{opts}, 
       FilterRules[Options[
         Graphics], Except[
         {ImagePadding, Frame, 
          FrameTicks}]]]}, 
    {ny, nx} = Dimensions[l]; 
     widths = 
      ((w - 2*sidePadding)/nx)*
       Table[1, {nx}]; 
     widths[[1]] = widths[[1]] + 
       sidePadding; 
     widths[[-1]] = 
      widths[[-1]] + sidePadding; 
     heights = 
      ((h - 2*sidePadding)/ny)*
       Table[1, {ny}]; 
     heights[[1]] = 
      heights[[1]] + sidePadding; 
     heights[[-1]] = 
      heights[[-1]] + 
       sidePadding; positions = 
      Transpose[Partition[
        Tuples[(Prepend[
           Accumulate[Most[#1]], 
           0] & ) /@ {widths, 
           heights}], ny]]; 
     Graphics[Table[Inset[
        Show[l[[ny - j + 1,i]], 
         ImagePadding -> 
          {{If[i == 1, 
           sidePadding, 0], 
           If[i == nx, 
           sidePadding, 0]}, 
           {If[j == 1, 
           sidePadding, 0], 
           If[j == ny, 
           sidePadding, 
           topPadding]}}, 
         AspectRatio -> Full], 
        positions[[j,i]], 
        {Left, Bottom}, 
        {widths[[i]], heights[[
          j]]}], {i, 1, nx}, 
       {j, 1, ny}], PlotRange -> 
       {{0, w}, {0, h}}, 
      ImageSize -> {w, h}, 
      Evaluate[Sequence@@frameOptions]]]; 


(* ::Input::Initialization:: *)
End[]


(* ::Input::Initialization:: *)
EndPackage[]
