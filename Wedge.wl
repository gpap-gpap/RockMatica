(* ::Package:: *)

(*  Mathematica Package  *)
(* :Title:   RockMatica`Wedge`  *)
(* :Author:  Giorgos Papageorgiou <gpap.gpap@gmail.com> *)
(* :Context: RockMatica`Wedge` *)
(* :Date:    01/11/2019  *)

(* :Version: 12 *)
(* :Copyright: (c) 2019 Giorgos Papageorgiou *)


(* ::Input::Initialization:: *)
BeginPackage["RockMatica`Wedge`"]


(* ::Input::Initialization:: *)
rpWedgeWavelet::usage="rpWedgeWavelet[f, f0] returns the spectral value of a ricker wavelet of central frequency f0, at frequency f";
rpWedgeSymConjAr::usage= "rpWedgeSymConjAr[FreqArray] symmetrizes and conjugates the array FreqArray in a form to be used by the inverse fourier transform. It assumes the first element corresponds to zero value if the array has odd length, otherwise adds  zero value";
rpWedgeInvFourier::usage="rpWedgeInvFourier[SymConjArray] inverse Fourier transforms the array SymConjArray";
rpWedgeElastWedgeRef::usage="rpWedgeElastWedgeRef[maxTime, res, waveletf, opts] creates an array with dimensions (2*res, 2*maxSamples + 1) corresponding to 1D convolution of a Ricker wavelet with a wedge model (first res lines correpond to even, last res lines to odd wedge)";
rpWedgeDispWedgeRef::usage="rpWedgeDispWedgeRef[rTop, rBot, a, maxTime, res, waveletf, opts] creates an array with dimensions (`res`, `2*maxSamples + 1`)  corresponding to 1D convolution of a Ricker wavelet with a dispersive wedge model with top reflectivity `rTop`, bottom reflectivity `rBot` and attenuation `(a-1)/(2 Sqrt[a])` (see Papageorgiou and Chapman; 2020).";
rpWedgeModelElasticTraces::usage = "rpWedgeModelElasticTraces[NReflections, nTraces, fwav, maxAbsRef, opts] creates an array with dimensions (`nTraces`, `2*maxSamples + 1`) corresponding to `nTraces` with at most `NReflections` each, convolved with a ricker wavelet of frequency `fwav`. The reflectivities and locations of traces are chosen at random which is repeatable by setting the random seed option `seed`";
rpWedgeModelReflectivities::usage = "rpWedgeModelReflectivities[NReflections, nTraces] creates an array with dimensions (`nTraces`, `2*maxSamples + 1`) corresponding to `nTraces` geological models with at most `NReflections` spikes each. This is essentially the deconvolved reflectivity series used to create elastic traces `rpWedgeCreateElasticTraces` if called with same number of traces and reflections.";
rpWedgeModelViscoelasticTraces::usage="rpWedgeModelViscoelasticTraces[sparseRef, maxThickness, waveletfreq, attenuation] forward models the layer model `sparseRef` with layers satisfying the thickness threshold `maxThickness` being viscoelastic. Due to limitations in the theoretical model that assumes elastic halfspaces enclosing the wedge (Papageorgiou and Chapman; 2020), the model looks at odd numbered layers unlesss the option `replaceOdd` is set to False.";
rpWedgeCreateNoise::usage="rpWedgeCreateNoise[nTraces, freq] creates noise bandlimited up to frequency freq (by default set to twice the wavelet frequency).";
rpInitiateThinLayerDictionary::usage="Create thin layer dictionary";
rpInitiateWaveletDictionary::usage="Create thin layer dictionary";
rpInitiatePhasedWaveletDictionary::usage="Create phased wavelet dictionary";
rpSparseVectorToSpikes::usage="";
rpWedgeBasisPursuitInversion::usage="";
rpWedgeDictionaryUpdate::usage="";
rpWedgeCompareResiduals::usage="";


(* ::Input::Initialization:: *)
Begin["Private`"]


(*default options*)
$RandomSeed = 1234567890;
$SampleRate = 0.001;
$MaxSamples = 2^11;
$scaleFactors = {1., 1.};
(*default values*)
$$WaveletFrequency=35;
$$MinQ = 1.22;
(*default functions*)
reflSeries[rTop_, rBot_, fr_, fwav_, alpha_, dt_]:= rTop + (1. - (fr/fwav)^2 - 2*I*fr/fwav )/(8.(1+(fr/fwav)^2)) (rTop^2-1)(alpha-1) + (rBot - (1. - (fr/fwav)^2 - 2*I*fr/fwav )/(8.(1+(fr/fwav)^2)) (rBot^2-1)(alpha-1)) E^(2.*I*Pi*fr*dt)//Chop//N;


rpWedgeWavelet[f_, f0:(_?NumericQ):$$WaveletFrequency] := With[{omega = 2.*Pi*f, omega0 = 2.*Pi*f0},
2.*Exp[-omega^2/omega0^2]*(omega^2/omega0^3)/Sqrt[Pi]];


rpWedgeSymConjAr[freqArray_?ArrayQ]:=With[{revconj=Reverse@Conjugate@freqArray},
Which[
	Mod[Length@freqArray,2]==1,Most@revconj~Join~freqArray,
	Mod[Length@freqArray,2]==0,revconj~Join~{0.}~Join~freqArray
]
];


rpWedgeInvFourier[freqArray_?ArrayQ]:=With[{n=Length@freqArray, Params={1,1}},
Which[
	Mod[n,2]==1, RotateRight[InverseFourier[RotateLeft[freqArray,(n-1)/2], FourierParameters->Params],(n-1)/2],
	Mod[n,2]==0, RotateLeft[InverseFourier[RotateRight[freqArray,(n)/2], FourierParameters->Params],(n)/2]
]
];


Options[rpWedgeElastWedgeRef] = {sampleRate -> $SampleRate, maxSamples -> $MaxSamples, scaleFactors->$scaleFactors, scaleType->"biased"}
rpWedgeElastWedgeRef[maxT:(_?NumericQ):1/$$WaveletFrequency, res:(_?NumericQ):20, fwav:(_?NumericQ):$$WaveletFrequency, OptionsPattern[]]:=
With[{fr=Subdivide[0., 1/(2.*OptionValue[sampleRate]), OptionValue[maxSamples]], invFour = rpWedgeInvFourier@rpWedgeSymConjAr@#&},
	Module[{reflSeriesOdd, reflSeriesEven, waveSingle, waveEven, waveOdd, resc = (#/Max@Abs@N@Flatten@#)&},
		reflSeriesOdd[f_, dt_] := 1. - E^(2.*I*Pi*f*dt);
		reflSeriesEven[f_, dt_] := 1. + E^(2.*I*Pi*f*dt);
			With[{wavelet = N@Chop@rpWedgeWavelet[fr,fwav], scFac = OptionValue[scaleFactors]},
				Which[OptionValue[scaleType]=="biased",
					waveSingle = scFac[[1]]*resc@N@Chop@invFour@(wavelet reflSeriesEven[fr,0.]);
					waveEven = scFac[[2]]*resc@Table[N@Chop@invFour@(wavelet*reflSeriesEven[fr,dt]),{dt, Rest@Subdivide[maxT/res, maxT, res-1]}];
					waveOdd = scFac[[2]]*resc@Table[N@Chop@invFour@(wavelet*reflSeriesOdd[fr,dt]), {dt, Rest@Subdivide[maxT, res]}];
					,
					OptionValue[scaleType]=="unbiased",
					waveSingle = scFac[[1]]*resc@N@Chop@invFour@(wavelet reflSeriesEven[fr,0.]);
					waveEven = Table[scFac[[2]]*resc@N@Chop@invFour@(wavelet*reflSeriesEven[fr,dt]),{dt, Rest@Subdivide[maxT/res, maxT, res-1]}];
					waveOdd = Table[scFac[[2]]*resc@N@Chop@invFour@(wavelet*reflSeriesOdd[fr,dt]), {dt, Rest@Subdivide[maxT, res]}];
					];
			Developer`ToPackedArray[{waveSingle}~Join~waveEven~Join~waveOdd]
		]
	]
];


Options[rpWedgeDispWedgeRef] = {sampleRate -> $SampleRate, maxSamples -> $MaxSamples, scaleType->"biased"}
rpWedgeDispWedgeRef[rTop:(_?NumericQ), rBot:(_?NumericQ), alpha:(_?NumericQ), maxT:(_?NumericQ):1/$$WaveletFrequency, res:(_?NumericQ):8, fwav:(_?NumericQ):$$WaveletFrequency, OptionsPattern[]]:=
With[{fr=Subdivide[0., 1/(2.*OptionValue[sampleRate]), OptionValue[maxSamples]], invFour = rpWedgeInvFourier@rpWedgeSymConjAr@#&},
	Module[{ref, waveForm, resc = (Max[Abs[rTop],Abs[rBot]]#/Max@Abs@Flatten@#)&},
		ref[dt_]:= reflSeries[rTop, rBot, fr, fwav, alpha, dt];
			With[{wavelet = N@Chop@rpWedgeWavelet[fr,fwav]},
				Which[OptionValue[scaleType]=="biased",
					waveForm = resc@Table[N@Chop@invFour@(wavelet*ref[dt]),{dt, Rest@Subdivide[maxT, res]}];
					,
					OptionValue[scaleType]=="unbiased",
					waveForm = Table[resc@N@Chop@invFour@(wavelet*ref[dt]),{dt, Rest@Subdivide[maxT, res]}];
					];
			Developer`ToPackedArray[waveForm]
		]
	]
];


Options[rpWedgeModelElasticTraces]={seed->$RandomSeed, sampleRate -> $SampleRate, maxSamples -> $MaxSamples, cutoffAmplitude->10^(-3)}
rpWedgeModelElasticTraces[NReflections_, nTraces_, fwav:(_?NumericQ):$$WaveletFrequency, maxAbsRef_Real:0.2, opts:OptionsPattern[]]:=
With[{fr=Subdivide[0., 1/(2.*OptionValue[sampleRate]), OptionValue[maxSamples]]},
	Module[{positions, magnitudes, Wlet, thread},
			SeedRandom[OptionValue[seed]];
			Wlet = (*(#/Max@Abs@#)&*)Normalize@rpWedgeInvFourier@rpWedgeSymConjAr@N@Chop@rpWedgeWavelet[fr,fwav];
			positions = RandomChoice[Range[0,OptionValue[maxSamples]], nTraces*NReflections];
			magnitudes = RandomReal[{-maxAbsRef, maxAbsRef}, nTraces*NReflections];
			thread = MapThread[#2 RotateRight[Wlet,#1]&,{positions, magnitudes}]//Chop//N;
			Plus@@@Partition[Threshold[thread, OptionValue[cutoffAmplitude]*maxAbsRef][[;;,OptionValue[maxSamples]+1;;]], NReflections]
		]
]


Options[rpWedgeModelReflectivities]={seed->$RandomSeed, sampleRate->$SampleRate, maxSamples->$MaxSamples}
rpWedgeModelReflectivities[NReflections_, nTraces_, maxAbsRef_Real:0.2, OptionsPattern[]]:=
Module[{positions, magnitudes, thread},
	SeedRandom[OptionValue[seed]];
	positions = RandomChoice[Range[0,OptionValue[maxSamples]], nTraces*NReflections];
	magnitudes = RandomReal[{-maxAbsRef, maxAbsRef}, nTraces*NReflections];
	thread = MapThread[SparseArray[(#1+1)->#2, OptionValue[maxSamples]+1]&,{positions, magnitudes}];
	Plus@@@Partition[thread, NReflections]
]


Options[rpWedgeModelViscoelasticTraces]={sampleRate -> $SampleRate, maxSamples -> $MaxSamples, replaceOdd->True};
rpWedgeModelViscoelasticTraces[a_SparseArray, crit:(_?NumericQ):1/$$WaveletFrequency, fwav:(_?NumericQ):$$WaveletFrequency, atten:(_?NumericQ):$$MinQ, opts:OptionsPattern[]]:=
With[{fr=Subdivide[0., 1/(2.*OptionValue[sampleRate]), OptionValue[maxSamples]],
invFour = rpWedgeInvFourier@rpWedgeSymConjAr@#&, 
sRate=OptionValue[sampleRate], rules=Most@ArrayRules[a], 
mxsamp=OptionValue[maxSamples], replOdd=OptionValue[replaceOdd]},
	Module[{disp, el, dispersify, intList},
	
	(*Define dispersive and elastic reflectivity for a layer*)
		disp[r1_, r2_, dt1_, dt2_]:=(*(Max[Abs[r1],Abs[r2]]#/Max@Abs@#)&@*)invFour[E^(I*2*Pi*fr*(dt1-1)*sRate)*reflSeries[r1, r2, fr, fwav, atten, (dt2-dt1)*sRate]*rpWedgeWavelet[fr,fwav]];
		el[r1_, r2_, dt1_, dt2_]:=  (*(Max[Abs[r1],Abs[r2]]#/Max@Abs@#)&@*)invFour[(r1*E^(I*2*Pi*fr*(dt1-1)*sRate)+r2*E^(I*2*Pi*fr*(dt2-1)*sRate))*rpWedgeWavelet[fr,fwav]];
	
	(*Define a function that chooses the dispersive reflectivity
	if the layer is thinner than 'crit'*)
		dispersify[{pos1_}->r1_, {pos2_}->r2_]:=Which[
				(pos2-pos1)*sRate <= crit, disp[r1, r2, pos1, pos2],
				(pos2-pos1)*sRate > crit, el[r1, r2, pos1, pos2]
			];
		dispersify[{pos1_}->r1_]:= (*(Abs[r1]#/Max@Abs@#)&@*)invFour[r1*E^(I*2*Pi*fr*(pos1-1)*sRate)*rpWedgeWavelet[fr,fwav]];
	
	(*Define a function that groups reflection pairs accounting for
	odd number of layers. Based on boolean option 'replOdd' make
	odd numbered (if replOdd is True) or even numbered (if replOdd 
	is False) layers potentially dispersive*)
		Which[
				Length@rules==1,intList = {rules} ,
				(Mod[Length@rules,2] == 0)&&replOdd, intList = Partition[rules,2],
				(Mod[Length@rules,2] == 0)&&(!replOdd), intList = {{First@rules}}~Join~Partition[Most@Rest@rules,2]~Join~{{Last@rules}},
				(Mod[Length@rules,2] == 1)&&replOdd, intList = Partition[Most@rules,2]~Join~{{Last@rules}},
				(Mod[Length@rules,2] == 1)&&(!replOdd), intList = {{First@rules}}~Join~Partition[Rest@rules,2]
			];
			
	(*Apply function that makes layers dispersive to every pair 
	of layers prepared in intList. Then add all the layers together
	to obtain the new trace*)
			(*Max@Abs@a(#/Max@Abs@#)&@*)(Plus@@(dispersify@@@intList))[[mxsamp + 1 ;;]]//Chop//N
	]
]


Options[rpWedgeCreateNoise]={seed->$RandomSeed, sampleRate -> $SampleRate, maxSamples -> $MaxSamples}
rpWedgeCreateNoise[nTraces_, fLimit:(_?NumericQ):2*$$WaveletFrequency, maxLevel_:0.01, OptionsPattern[]]:=
	Module[{fourierSeries, bandLimitNoise},
			SeedRandom[OptionValue[seed]];
			fourierSeries = RandomComplex[{-1-2*Pi*I, 1+2*Pi*I}, Ceiling[2*Pi*fLimit]];
			bandLimitNoise = PadRight[fourierSeries, OptionValue[maxSamples]];
			(maxLevel #/Max@Abs@#)&@(rpWedgeInvFourier@rpWedgeSymConjAr[bandLimitNoise])[[OptionValue[maxSamples] +  1 ;;]]//Chop//N
		]


rpWedgeBasisPursuitInversion[dictionary_Array, dataVector:{_?NumericQ..}, regParam_Real/;(regParam > 0), opts:OptionsPattern[]]:=
Which[
	First@Dimensions@dictionary == Length@dataVector, Fit[{dictionary, dataVector}, FilterRules[opts,Options[Fit]], FitRegularization->{"L1",regParam}],
	True, $Failed
]



Options[rpInitiateWaveletDictionary]={sampleRate->$SampleRate, maxSamples->$MaxSamples, threshold->0.001, scaleType->"unbiased"};
rpInitiateWaveletDictionary[freq:(_?NumericQ):$$WaveletFrequency, opts:OptionsPattern[]]:=
With[{fr=Subdivide[0., 1/(2.*OptionValue[sampleRate]), OptionValue[maxSamples]], invFour = rpWedgeInvFourier@rpWedgeSymConjAr@#&,thres=OptionValue[threshold]},
	With[{wlet=(*(#/Max@Abs@N@#)&*)Normalize@Chop@invFour@rpWedgeWavelet[fr,freq]},
		With[{rot=(RotateRight[wlet, #][[OptionValue[maxSamples]+1;;]])&/@Range[0,OptionValue[maxSamples]]},
				SparseArray@Threshold[rot,thres]
		]
	]
]


Options[rpInitiateThinLayerDictionary]={sampleRate->$SampleRate, maxSamples->$MaxSamples, threshold->0.001, scaleType->"unbiased"};
rpInitiateThinLayerDictionary[maxT_, tSteps_Integer, freq:(_?NumericQ):$$WaveletFrequency, sSteps_Integer:1, opts:OptionsPattern[]]:=
With[{w=rpWedgeElastWedgeRef[maxT, tSteps, freq, FilterRules[{opts},Options[rpWedgeElastWedgeRef]]], thres=OptionValue[threshold], samp=OptionValue[maxSamples]},
	With[{wEven = Normalize/@w[[ ;; tSteps]], wOdd = Normalize/@w[[tSteps + 1 ;;]]},
		With[{rotEven=Table[RotateRight[i,j][[samp+1;;]],{j,0, samp, sSteps},{i,wEven}],
			   rotOdd=Table[RotateRight[i,j][[samp+1;;]],{j,0, samp, sSteps},{i,wOdd}]},
			SparseArray@Transpose@Threshold[Flatten[rotEven,1]~Join~Flatten[rotOdd,1],thres]
		]
	]	
]


Options[rpInitiatePhasedWaveletDictionary]={sampleRate->$SampleRate, maxSamples->$MaxSamples, threshold->0.0001, scaleType->"unbiased"};
rpInitiatePhasedWaveletDictionary[maxAbsFreq:(_?NumericQ), freq_:$$WaveletFrequency, opts:OptionsPattern[]]:=
With[{sRate=OptionValue[sampleRate], nSamps=OptionValue[maxSamples],thres=OptionValue[threshold]},
With[{fr=Subdivide[0.,1/(2 sRate), nSamps]},
	SparseArray@Transpose@
	Threshold[
	Flatten[
	Table[
		With[{wlet=Normalize@Chop@rpWedgeInvFourier@rpWedgeSymConjAr@(E^(I \[Phi])*rpWedgeWavelet[fr,freq])},
			With[{rot=(RotateRight[wlet, #][[nSamps+1;;]])&/@Range[0,nSamps]},
				rot
				]
			]
		,{\[Phi],{-maxAbsFreq,0, maxAbsFreq}}]
		,1]
		,thres]
	]
]


Options[rpSparseVectorToSpikes]={sampleRate->$SampleRate, maxSamples->$MaxSamples, scaleFactors->$scaleFactors};
rpSparseVectorToSpikes[sparsevector_SparseArray, maxT_, tSteps_Integer,opts:OptionsPattern[]]:=
With[{samps=OptionValue[maxSamples], sRate=OptionValue[sampleRate],scales=OptionValue[scaleFactors]},	
	With[{
		even = sparsevector[[;;(samps+1)*tSteps]], 
		odd =  sparsevector[[(samps+1)*tSteps+1;;]]},
			scales[[1]]*even[[1;; ;;tSteps]]+
			scales[[2]]*Plus@@Table[
				even[[i+1;; ;;tSteps]] + 
				RotateRight[even[[i+1;; ;;tSteps]], Floor[i*maxT/tSteps/sRate]],
			{i,1,tSteps-1}]+
			scales[[2]]*Plus@@Table[
				odd[[i;;;;tSteps]] - 
				RotateRight[odd[[i;; ;;tSteps]], Floor[i*maxT/tSteps/sRate]],
			{i,1,tSteps}]
		
		
	]
]


Options[rpWedgeDictionaryUpdate]={sampleRate -> $SampleRate}
rpWedgeDictionaryUpdate[dictionary_Array, dataVector:{_?NumericQ..}]:=Null;


rpWedgeCompareResiduals[result_/;ArrayQ@result, {vector1_/;ArrayQ@vector1, vector2_/;ArrayQ@vector2}]:=With[{f=Tr@Abs@#&},
		{f[result-vector1], f[result-vector2]}
	];


(* ::Input::Initialization:: *)
End[]


(* ::Input::Initialization:: *)
EndPackage[]
