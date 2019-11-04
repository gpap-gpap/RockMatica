(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["RockMatica`Moduli`"]


(* ::Input::Initialization:: *)
rpIsotropicElasticTensor::usage="rpIsotropicElasticTensor[a_rpRock, b_rpFluid] returns the elastic tensor (in Voigt notation) of a rock based on Gassmann's formula. If additionally the option 'Viscoelastic->True, Frequency->freq' are given, a viscoelastic, complex-valued tensor is returned based on the Chapman et al. (2002) model. In that case, the rpRock object 'a' needs to have defined crack density, aspect ratio and reference frequency predefined.";
rpTIAnisotropicElasticTensor::usage = "rpTIAnisotropicElasticTensor[a_rpRock, b_rpFluid] returns the VTI anisotropic elastic tensor (in Voigt notation) of a rock based on the anisotropic Gassmann's formula. If additionally the option 'Viscoelastic->True, Frequency->freq' are given, a viscoelastic, complex-valued tensor is returned based on the Chapman et al. (2002) model. In that case, the rpRock object 'a'needs to have defined crack density, fracture density, aspect ratio, reference frequency and ratio of fracture-to-microcrack reference frequency predefined.";


(* ::Input::Initialization:: *)
Begin["Private`"]


(* Gassmann's isotropic and TI anisotropic elastic tensors *)
gassmannIsotropic[Kd_, mu_, Km_, phi_, Kf_]:=Module[{Kgas},
	Kgas = Kd+(1-Kd/Km)^2/(-(Kd/Km^2)+(1-phi)/Km+phi/Kf);
	SparseArray[{
	{1,1} -> Kgas+4./3.mu,
	{2,2} -> Kgas+4./3.mu,
	{3,3} -> Kgas+4./3.mu,
	{1,2} -> Kgas-2./3.mu,
	{2,1} -> Kgas-2./3.mu, 
	{1,3} -> Kgas-2./3.mu,
	{2,3} -> Kgas-2./3.mu,
	{3,2} -> Kgas-2./3.mu,
	{3,1} -> Kgas-2./3.mu,
	{4,4} -> mu,
	{5,5} -> mu,
	{6,6} -> mu}]
    ];
gassmannTIAnisotropic[c11d_, c12d_, c13d_, c33d_, c44d_, Km_, phi_, Kf_]:=Module[{norm},
norm = Km/Kf phi (Km - Kf) + Km -(2 c11d + 2 c12d + 4c13d + c33d)/9.;
SparseArray[{
	{1,1} -> c11d + (Km - (c11d +c12d +c13d)/3)^2/norm,
	{2,2} -> c11d + (Km - (c11d +c12d +c13d)/3)^2/norm,
	{3,3} -> c33d + (Km - (c33d +c13d +c13d)/3)^2/norm,
	{1,2} -> c12d + (Km - (c11d +c12d +c13d)/3)^2/norm,
	{2,1} -> c12d + (Km - (c11d +c12d +c13d)/3)^2/norm, 
	{1,3} -> c13d + (Km - (c11d +c12d +c13d)/3)(Km - (c33d +c13d +c13d)/3)/norm,
	{2,3} -> c13d + (Km - (c11d +c12d +c13d)/3)(Km - (c33d +c13d +c13d)/3)/norm,
	{3,2} -> c13d + (Km - (c11d +c12d +c13d)/3)(Km - (c33d +c13d +c13d)/3)/norm,
	{3,1} -> c13d + (Km - (c11d +c12d +c13d)/3)(Km - (c33d +c13d +c13d)/3)/norm,
	{4,4} -> c44d,
	{5,5} -> c44d,
	{6,6} -> c11d + (3Km - (c11d +c12d +c13d)/3)^2/norm-c12d + (3Km - (c11d +c12d +c13d)/3)^2/norm}]
];


(* Squirt isotropic and TI anisotropic elastic tensors *)
squirtIsotropic[mud_, Kd_, Km_, phi_, Kf_, epsilon_, aspectRatio_, tau_, omega_] := Module[
    {Kc, Kp, gamma, gamma2, nu, sigmac, phic, phip, muhf, mu, Khf, Kgassmann, Ksquirt, musquirt},
    phic = 4./3. Pi epsilon aspectRatio; 
	phip = phi - phic;
	nu =(3.Km - 2.mu)/(2(3Km + mu));
	mu =-((16. Km^2 phic+3 aspectRatio Pi Km (4 Kd+Km (-4+5 phip))+\[Sqrt](Km^2 (256 Km^2 phic^2+9 aspectRatio^2 Pi^2 (4 Kd-4 Km+3 Km phip)^2+96 aspectRatio Pi Km phic (2 Kd-2 Km+3 Km phip))))/(8 aspectRatio Pi (Kd+Km (-1+phip))));
	sigmac = (3. aspectRatio \[Pi] Km (-1+2 nu))/(4 (-1+nu^2)); 
	Kc = sigmac/Kf; 
	Kp =(4.*mu)/(3*Kf); 
	gamma =(3.*phip*sigmac*(1 + Kp))/(4*phic*mu*(1 + Kc)); 
    gamma2 = gamma*((1 - nu)/((1 + nu)*(1 + Kp))); 
	Kgassmann = Kd+(Km (1+3 (1+Kc) gamma2) ((1+Km/sigmac) phic+(1+(3 Km)/(4 mu)) phip))/((1+Kc) (1+gamma)); 
	Khf =(Km (gamma-3 (1+Kc) gamma2) (gamma (1+Km/sigmac) phic-(1+(3 Km)/(4 mu)) phip))/((1+Kc) gamma (1+gamma) (1-(I (1+gamma))/(gamma tau omega)));
	muhf = (4 mu^2 omega phic tau)/(15 (1+Kc) sigmac (-I+omega tau))(* *);
	Ksquirt = Kgassmann+Khf;
	musquirt = mud + muhf;
	SparseArray[{
	{1,1} -> Ksquirt+4./3.musquirt,
	{2,2} -> Ksquirt+4./3.musquirt,
	{3,3} -> Ksquirt+4./3.musquirt,
	{1,2} -> Ksquirt-2./3.musquirt,
	{2,1} -> Ksquirt-2./3.musquirt, 
	{1,3} -> Ksquirt-2./3.musquirt,
	{2,3} -> Ksquirt-2./3.musquirt,
	{3,2} -> Ksquirt-2./3.musquirt,
	{3,1} -> Ksquirt-2./3.musquirt,
	{4,4} -> musquirt,
	{5,5} -> musquirt,
	{6,6} -> musquirt}]
    ]
    
squirtTIAnisotropic[lambda_, mu_, phi_, Kf_, epsilon_, epsilonf_, aspectRatio_, tau0_, omega_, lenrat_]:=Module[
	{phic,phif,phip,gamma,gamma2,iota,beta,sigmac,nu,Kp,Kc,D1,D2,G1,G2,G3,F1,F2,k,L2,L4,tauf,taum,c},
	phic=(4. \[Pi] epsilon aspectRatio)/3.;
	phif=(4. \[Pi] epsilonf aspectRatio)/3.;
	phip=phi-phic-phif;
	Kp=(4. mu)/(3. Kf);
	Kc=sigmac/Kf;
	nu=lambda/(2. (lambda+mu));
	sigmac=(\[Pi] mu aspectRatio)/(2. (1-nu));
	phip=phi-phic-phif;
	gamma=(3. phip sigmac (1+Kp))/(4. phic mu (1+Kc));
	gamma2=(gamma (1.-nu))/((1.+nu) (1.+Kp));
	iota=phic/(phic+aspectRatio phip);
	beta=(iota phif)/phic;
	tauf=lenrat tau0;
	taum=tau0;
	D1=(iota (-I+tauf omega-beta taum omega)+3 (1+Kc) gamma2 ((1+I tauf omega) (-I+taum omega)+iota (I-tauf omega+beta taum omega)))/(3 (1+Kc) (beta (-I+(1+(-1+gamma) iota) taum omega)-(-I+tauf omega) (-iota+gamma (-1+iota-I taum omega))));
	D2=(beta (-I+taum omega))/((1+Kc) (beta (-I+(1+(-1+gamma) iota) taum omega)-(-I+tauf omega) (-iota+gamma (-1+iota-I taum omega))));
	G1=(taum omega)/((1+Kc) (-I+taum omega));
	G2=(I (I-tauf omega+beta taum omega) (iota+I gamma iota taum omega-3 I (1+Kc) gamma2 (-1+iota) (-I+taum omega)))/(3 (1+Kc) (-I+taum omega) (beta (-I+(1+(-1+gamma) iota) taum omega)-(-I+tauf omega) (-iota+gamma (-1+iota-I taum omega))));
	G3=(beta (-I+gamma taum omega))/((1+Kc) (beta (-I+(1+(-1+gamma) iota) taum omega)-(-I+tauf omega) (-iota+gamma (-1+iota-I taum omega))));
	F1=(-3 (1+Kc) gamma2 (-1+iota) (-I+taum omega)+iota (-I+gamma taum omega))/(3 (1+Kc) (beta (-I+(1+(-1+gamma) iota) taum omega)-(-I+tauf omega) (-iota+gamma (-1+iota-I taum omega))));
	F2=(tauf omega (gamma+iota-gamma iota+I gamma taum omega)+beta (-I+(1+(-1+gamma) iota) taum omega))/((1+Kc) (beta (-I+(1+(-1+gamma) iota) taum omega)-(-I+tauf omega) (-iota+gamma (-1+iota-I taum omega))));
	k=lambda+(2 mu)/3;
	L2=k^2+(16 mu^2)/45;
	L4=k^2-(8 mu^2)/45;
	c[11]=lambda+2 mu-
		phic (L2/sigmac+(32 (1-nu) mu)/(15 (2-nu) (\[Pi] aspectRatio))-(L2/sigmac+k) G1-((3 k^2)/sigmac+3 k) G2-((lambda k)/sigmac+lambda) G3)-
		phip ((3 (1-nu) (3 lambda^2+4 lambda mu+((36+20 nu) mu^2)/(7-5 nu)))/((4 mu) (1+nu))-(1+(3 k)/(4 mu)) (3 k D1+lambda D2))-
		phif (lambda^2/sigmac-3 k (lambda/sigmac+1) F1-lambda (lambda/sigmac+1) F2);
	c[33]=lambda+2 mu-
		phic (L2/sigmac+(32 (1-nu) mu)/(15 (2-nu) (\[Pi] aspectRatio))-(L2/sigmac+k) G1-((3 k^2)/sigmac+3 k) G2-(((lambda+2 mu) k)/sigmac+lambda+2 mu) G3)-
		phip ((3 (1-nu) (3 lambda^2+4 lambda mu+((36+20 nu) mu^2)/(7-5 nu)))/((4 mu) (1+nu))-(1+(3 k)/(4 mu)) (3 k D1+(lambda+2 mu) D2))-
		phif ((lambda+2 mu)^2/sigmac-3 k ((lambda+2 mu)/sigmac+1) F1-(lambda+2 mu) ((lambda+2 mu)/sigmac+1) F2);
	c[44]=mu-
		phic ((4 mu^2 (1-G1))/(15 sigmac)+(8 (1-nu) mu)/(5 (2-nu) (\[Pi] aspectRatio)))-
		(phip (15 mu (1-nu)))/(7-5 nu)-
		(phif (4 (1-nu) mu))/((2-nu) (\[Pi] aspectRatio));
	c[12]=lambda-
		phic (L4/sigmac-(16 (1-nu) mu)/(15 (2-nu) (\[Pi] aspectRatio))-(L4/sigmac+k) G1-((3 k^2)/sigmac+3 k) G2-((lambda k)/sigmac+lambda) G3)-
		phip ((3 (1-nu) (3 lambda^2+4 lambda mu-((4 (1+5 nu)) mu^2)/(7-5 nu)))/((4 mu) (1+nu))-(1+(3 k)/(4 mu)) (3 k D1+lambda D2))-
		phif (lambda^2/sigmac-3 k (lambda/sigmac+1) F1-lambda (lambda/sigmac+1) F2);
	c[13]=lambda-
		phic (L4/sigmac-(16 (1-nu) mu)/(15 (2-nu) (\[Pi] aspectRatio))-(L4/sigmac+k) G1-((3 k^2)/sigmac+3 k) G2-(((lambda+mu) k)/sigmac+lambda+mu) G3)-
		phip ((3 (1-nu) (3 lambda^2+4 lambda mu-((4 (1+5 nu)) mu^2)/(7-5 nu)))/((4 mu) (1+nu))-(1+(3 k)/(4 mu)) (3 k D1+(lambda+mu) D2))-
		phif ((lambda (lambda+2 mu))/sigmac-3 k ((lambda+mu)/sigmac+1) F1-((lambda (lambda+2 mu))/sigmac+lambda+mu) F2);
SparseArray[{
	{1,1} -> c[11],
	{2,2} -> c[11],
	{3,3} -> c[33],
	{1,2} -> c[12],
	{2,1} -> c[21], 
	{1,3} -> c[13],
	{2,3} -> c[13],
	{3,2} -> c[13], 
	{3,1} -> c[13],
	{4,4} -> c[44],
	{5,5} -> c[44],
	{6,6} -> c[11]-c[12]}]
];


Options[rpIsotropicElasticTensor] = {Viscoelastic->False, Frequency->None};
SetAttributes[rpIsotropicElasticTensor,HoldAll];
rpIsotropicElasticTensor[a_RockMatica`Base`rpRock, b_RockMatica`Base`rpFluid, OptionsPattern[]]:=
Module[{
	elasticIsoQ = Query[{#DryModulus, #ShearModulus, #MineralModulus, #Porosity}&],
	elasticTIQ = Query[{#C11, #C12, #C13, #C33, #C44, #MineralModulus, #Porosity}&],
	microcrackQ = Query[{#CrackDensity, #AspectRatio, #ReferenceFrequency}&],
	frackQ = Query[{#FractureDensity, #FractureFrequencyRatio}&],
	fluidQ = Query[#FluidModulus&]},
Check[
	Which[
	!OptionValue[Viscoelastic],
		Check[Which[KeyExistsQ["DryModulus"][a],
			Module[{Kd, Km, mu, phi,Kf},
				{Kd, mu, Km, phi}=elasticIsoQ[a];
				Kf=fluidQ[b];
				gassmannIsotropic[Kd, mu, Km, phi, Kf]
			],
			KeyExistsQ["C11"][a],
			Module[{C11, C12, C13, C33, C44, Km, phi,Kf},
				{C11, C12, C13, C33, C44, Km, phi}=elasticTIQ[a];
				Kf=fluidQ[b];
				gassmannTIAnisotropic[C11, C12, C13, C33, C44, Km, phi, Kf]
				]
			],
			$Failed
		],
	OptionValue[Viscoelastic], 
		Check[Which[
			KeyExistsQ["DryModulus"][a]&&(And@@(NumericQ/@microcrackQ[a["MicrocrackParameters"]])&&NumericQ@OptionValue[Frequency]),
			Module[{Kd, Km, mu,emc, aspRat, taufreq, phi,Kf},
				{Kd, mu, Km, phi}=elasticIsoQ[a];
				Kf=fluidQ[b];
				{emc, aspRat, taufreq}=microcrackQ[a["MicrocrackParameters"]];
				squirtIsotropic[mu, Kd, Km, phi,Kf, emc, aspRat, taufreq, OptionValue[Frequency]]
			],
			KeyExistsQ["DryModulus"][a]&&
			(And@@(NumericQ/@microcrackQ[a["MicrocrackParameters"]]))&&
			(And@@(NumericQ/@microcrackQ[a["FractureParameters"]])&&
			NumericQ@OptionValue[Frequency]),
			Module[{Kd, Km, mu,emc, aspRat, taufreq, phi,Kf},
				{Kd, mu, Km, phi}=elasticIsoQ[a];
				Kf=fluidQ[b];
				{emc, aspRat, taufreq}=microcrackQ[a["MicrocrackParameters"]];
				squirtIsotropic[mu, Kd, Km, phi,Kf, emc, aspRat, taufreq, OptionValue[Frequency]]
			],
			KeyExistsQ["C11"][a]&&
			(And@@(NumericQ/@microcrackQ[a["MicrocrackParameters"]]))&&
			(And@@(NumericQ/@microcrackQ[a["FractureParameters"]])&&
			NumericQ@OptionValue[Frequency]),
			Module[{C11, C12, C13, C33, C44, Km, phi,Kf},
				{C11, C12, C13, C33, C44, Km, phi}=elasticTIQ[a];
				Kf=First@fluidQ[b];
				{emc, aspRat, taufreq}=microcrackQ[a["MicrocrackParameters"]];
				{emf, len}=frackQ[a["FractureParameters"]];
					
				squirtTIAnisotropic[lam, mu, phi, Kf, emc, emf, aspRat, taufreq, OptionValue[Frequency], len]
				]
			],
			$Failed
		]
	],
		(*Module[{Kd, Km, mu, phi, Kf, etaf, emc, aspRat, taufreq},
			{Kd, mu, Km, phi}=elasticQ[a];
			{Kf, etaf}=fluidQ[b];
			{emc, aspRat, taufreq}=microcrackQ[a["MicrocrackParameters"]];
			squirtIsotropic[mu, Kd, Km, phi,Kf, emc, aspRat, taufreq, OptionValue[Frequency], etaf]
		]*)(*,
	(OptionValue[Viscoelastic]&&(And@@(NumericQ/@microcrackQ[a["MicrocrackParameters"]]))&&(And@@(NumericQ/@fractureQ[a["FractureParameters"]]))),
		Module[{Kd, Km, phi, Kf, etaf, emc, aspRat, taufreq, emf, relFreq},
			{Kd, Km, phi}=elasticQ[a];
			{Kf, etaf}=fluidQ[b];
			{emc, aspRat, taufreq}=microcrackQ[a["MicrocrackParameters"]];
			squirt2[Kd, Km, phi,Kf, emc, aspRat, taufreq, freq, etaf]
		]*)
$Failed]
]


(* ::Input::Initialization:: *)
End[]


(* ::Input::Initialization:: *)
EndPackage[]
