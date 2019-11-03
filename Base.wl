(* ::Package:: *)

BeginPackage["RockPhysics`Base`"]


rpRock::usage = "rpRock[RockName, DryBulkModulus, \
ShearModulus, MineralModulus, Porosity, MineralDensity]\
defines a rock object called by rpRock[RockName]. \
For an isotropically viscoelastic rock the option \
microcrackParameters->{CrackDensity, AspectRatio, ReferenceFrequency}\
can be given, whereas for an anisotropically viscoelastic rock\
fractureParameters->{FractureDensity, FractureFrequencyRatio}\
must also be provided
"


rpFluid::usage = "rpFluid[FluidName, FluidModulus, FluidViscosity, \
FluidDensity] defines a fluid object called by rpFluid[FluidName]"


rpAnisotropicRock::usage = 


Begin["Private`"]


Options[rpRock]={MicrocrackParameters->{None, 10^-5, None}, FractureParameters->{None,None}};
rpRock[name_String, drymod_, shearmod_, mineralmod_, porosity_, mineraldens_ , OptionsPattern[]]:=
Module[{crackDensity, aspectRatio, referenceFrequency, fractureDensity, fractureRefFreq},
{crackDensity, aspectRatio, referenceFrequency}=Evaluate@OptionValue[MicrocrackParameters];
{fractureDensity, fractureRefFreq}=Evaluate@OptionValue[FractureParameters];
rpRock[name]=<|
"RockName"->name, 
"DryModulus"->drymod, 
"ShearModulus"->shearmod,
"MineralModulus"->mineralmod, 
"Porosity"->porosity, 
"MineralDensity"->mineraldens, 
"MicrocrackParameters"-><|
	"CrackDensity"->crackDensity, 
	"AspectRatio"->aspectRatio,
	"ReferenceFrequency"->referenceFrequency
	|>,
"FractureParameters"-><|
	"FractureDensity"->fractureDensity,
	"FractureFrequencyRatio"->fractureRefFreq
	|>
|>//N
];
rpRock[name_String, c11_, c12_, c13_, c33_, c44_, mineralmod_, porosity_, mineraldens_ , OptionsPattern[]]:=
Module[{crackDensity, aspectRatio, referenceFrequency, fractureDensity, fractureRefFreq},
{crackDensity, aspectRatio, referenceFrequency}=Evaluate@OptionValue[MicrocrackParameters];
{fractureDensity, fractureRefFreq}=Evaluate@OptionValue[FractureParameters];
rpRock[name]=<|
"RockName"->name, 
"C11"->c11, 
"C12"->c12,
"C13"->c13,
"C33"->c33,
"C44"->c44,
"MineralModulus"->mineralmod, 
"Porosity"->porosity, 
"MineralDensity"->mineraldens, 
"MicrocrackParameters"-><|
	"CrackDensity"->crackDensity, 
	"AspectRatio"->aspectRatio,
	"ReferenceFrequency"->referenceFrequency
	|>,
"FractureParameters"-><|
	"FractureDensity"->fractureDensity,
	"FractureFrequencyRatio"->fractureRefFreq
	|>
|>//N
];


rpFluid[name_String,fluidMod_, fluidVisc_, fluidDens_]:=
rpFluid[name]=<|
"FluidName"->name, 
"FluidModulus"->fluidMod, 
"FluidViscosity"->fluidVisc,
"FluidDensity"->fluidDens
|>//N;


End[]


EndPackage[]
