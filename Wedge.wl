(* ::Package:: *)

(*  Mathematica Package  *)
(* :Title:   RockMatica  *)
(* :Author:  Giorgos Papageorgiou <gpap.gpap@gmail.com> *)
(* :Context: RockMatica`Wedge` *)
(* :Date:    01/11/2019  *)

(* :Version: 12 *)
(* :Copyright: (c) 2019 Giorgos Papageorgiou *)


(* ::Input::Initialization:: *)
BeginPackage["RockMatica`Wedge`"]


(* ::Input::Initialization:: *)
rpWedgeWavelet::usage="Create a ricker wavelet";
rpWedgeSymConjAr::usage= "Prepare array of values for inverse fourier transform (symmetrize, add zero value, conjugate";
rpWedgeInvFourier::usage="Inverse Fourier transform, taking into account the zero value";
rpWedgeElastWedgeRef::usage="Create an elastic wedge model";
rpWedgeDispWedgeRef::usage="Create a dispersive wedge model";
rpWedgeCreateElasticTraces::usage = "Create a random trace set";
rpWedgeCreateDispersiveTraces::usage = "Create a random trace set";
rpWedgeTracesReflections::usage = "Get reflections";
rpWedgeCreateNoise::usage="Create bandlimited noise";
rpThinLayerQ::usage="Detect Thin Layer in traces";
rpInitiateThinLayerDictionary::usage="Create thin layer dictionary";
(*rpWedgeFreqArray::usage=;
rpWedgeTimeArray::usage=;
*)


(* ::Input::Initialization:: *)
Begin["Private`"]


$RandomSeed = 1234567890;
$SampleRate = 0.001;
$MaxSamples = 2^11;


rpWedgeWavelet[f_, f0_] := With[{omega = 2.*Pi*f, omega0 = 2.*Pi*f0},
2.*Exp[-omega^2/omega0^2]*(omega^2/omega0^3)/Sqrt[Pi]];


rpWedgeSymConjAr[freqArray_]:=With[{revconj=Reverse@Conjugate@freqArray},
Which[
	Mod[Length@freqArray,2]==1,Most@revconj~Join~freqArray,
	Mod[Length@freqArray,2]==0,revconj~Join~{0.}~Join~freqArray
]
];


rpWedgeInvFourier[freqArray_]:=With[{n=Length@freqArray, Params={1,1}},
Which[
	Mod[n,2]==1, RotateRight[InverseFourier[RotateLeft[freqArray,(n-1)/2], FourierParameters->Params],(n-1)/2],
	Mod[n,2]==0, RotateLeft[InverseFourier[RotateRight[freqArray,(n)/2], FourierParameters->Params],(n)/2]
]
];


Options[rpWedgeElastWedgeRef] = {sampleRate -> $SampleRate, maxSamples -> $MaxSamples, scaleFactors->{1., 0.9}, scaleType->"biased"}
rpWedgeElastWedgeRef[maxT_, res_, fwav_, OptionsPattern[]]:=
With[{fr=Subdivide[0., 1/(2.*OptionValue[sampleRate]), OptionValue[maxSamples]], invFour = rpWedgeInvFourier@rpWedgeSymConjAr@#&},
	Module[{reflSeriesOdd, reflSeriesEven, waveSingle, waveEven, waveOdd, resc = (#/Max@Abs@N@Flatten@#)&},
		reflSeriesOdd[f_, dt_] := 1. - E^(2.*I*Pi*f*dt);
		reflSeriesEven[f_, dt_] := 1. + E^(2.*I*Pi*f*dt);
			With[{wavelet = N@Chop@rpWedgeWavelet[fr,fwav], scFac = OptionValue[scaleFactors]},
				Which[OptionValue[scaleType]=="biased",
					waveSingle = scFac[[1]]*resc@N@Chop@invFour@(wavelet reflSeriesEven[fr,0.]);
					waveEven = scFac[[2]]*resc@Table[N@Chop@invFour@(wavelet*reflSeriesEven[fr,dt]),{dt, Rest@Subdivide[maxT, res]}];
					waveOdd = scFac[[2]]*resc@Table[N@Chop@invFour@(wavelet*reflSeriesOdd[fr,dt]), {dt, Rest@Subdivide[maxT, res]}];
					,
					OptionValue[scaleType]=="unbiased",
					waveSingle = scFac[[1]]*resc@N@Chop@invFour@(wavelet reflSeriesEven[fr,0.]);
					waveEven = Table[scFac[[2]]*N@Chop@resc@invFour@(wavelet*reflSeriesEven[fr,dt]),{dt, Rest@Subdivide[maxT, res]}];
					waveOdd = Table[scFac[[2]]*N@Chop@resc@invFour@(wavelet*reflSeriesOdd[fr,dt]), {dt, Rest@Subdivide[maxT, res]}];
					];
			Developer`ToPackedArray[{waveSingle}~Join~waveEven~Join~waveOdd]
		]
	]
];


Options[rpWedgeDispWedgeRef] = {sampleRate -> $SampleRate, maxSamples -> $MaxSamples, scaleType->"biased"}
rpWedgeDispWedgeRef[rTop_, rBot_, alpha_, maxT_, res_, fwav_, OptionsPattern[]]:=
With[{fr=Subdivide[0., 1/(2.*OptionValue[sampleRate]), OptionValue[maxSamples]], invFour = rpWedgeInvFourier@rpWedgeSymConjAr@#&},
	Module[{reflSeries, waveForm, resc = (Max[Abs[rTop],Abs[rBot]]#/Max@Abs@Flatten@#)&},
		reflSeries[dt_]:= rTop + (1. - (fr/fwav)^2 - 2*I*fr/fwav )/(8.(1+(fr/fwav)^2)) (rTop^2-1)(alpha-1) + (rBot - (1. - (fr/fwav)^2 - 2*I*fr/fwav )/(8.(1+(fr/fwav)^2)) (rBot^2-1)(alpha-1)) E^(2.*I*Pi*fr*dt)//Chop//N;
			With[{wavelet = N@Chop@rpWedgeWavelet[fr,fwav]},
				Which[OptionValue[scaleType]=="biased",
					waveForm = resc@Table[N@Chop@invFour@(wavelet*reflSeries[dt]),{dt, Rest@Subdivide[maxT, res]}];
					,
					OptionValue[scaleType]=="unbiased",
					waveForm = Table[resc@N@Chop@invFour@(wavelet*reflSeries[dt]),{dt, Rest@Subdivide[maxT, res]}];
					];
			Developer`ToPackedArray[waveForm]
		]
	]
];


Options[rpWedgeCreateElasticTraces]={seed->$RandomSeed, sampleRate -> $SampleRate, maxSamples -> $MaxSamples}
rpWedgeCreateElasticTraces[NReflections_, nTraces_, fwav_, maxAbsRef_Real:0.2, OptionsPattern[]]:=
With[{fr=Subdivide[0., 1/(2.*OptionValue[sampleRate]), OptionValue[maxSamples]]},
	SeedRandom[OptionValue[seed]];
	(*Module[{positions, magnitudes, phasedWlet},
			phasedWlet = Function[{pos, amp}, amp(#/Max@Abs@#)&@rpWedgeInvFourier@rpWedgeSymConjAr@N@Chop@(Exp@(2*I*Pi*fr*pos)rpWedgeWavelet[fr,fwav])];
			positions = RandomChoice[Subdivide[OptionValue[sampleRate]*OptionValue[maxSamples],OptionValue[maxSamples]], nTraces*NReflections];
			magnitudes = RandomReal[{-maxAbsRef, maxAbsRef}, nTraces*NReflections];
			Plus@@@Partition[Threshold[MapThread[phasedWlet,{positions, magnitudes}]//Chop//N, 0.00001][[;;,OptionValue[maxSamples]+1;;]], NReflections]
		]*)
	Module[{positions, magnitudes, Wlet, thread},
			Wlet = (#/Max@Abs@#)&@rpWedgeInvFourier@rpWedgeSymConjAr@N@Chop@rpWedgeWavelet[fr,fwav];
			positions = RandomChoice[Range[0,OptionValue[maxSamples]], nTraces*NReflections];
			magnitudes = RandomReal[{-maxAbsRef, maxAbsRef}, nTraces*NReflections];
			thread = MapThread[#2 RotateRight[Wlet,#1]&,{positions, magnitudes}]//Chop//N;
			Plus@@@Partition[Threshold[thread, 0.00001][[;;,OptionValue[maxSamples]+1;;]], NReflections]
		]
]


Options[rpWedgeTracesReflections]={seed->$RandomSeed, sampleRate->$SampleRate, maxSamples->$MaxSamples}
rpWedgeTracesReflections[NReflections_, nTraces_, maxAbsRef_Real:0.2, OptionsPattern[]]:=
Module[{positions, magnitudes, thread},
	SeedRandom[OptionValue[seed]];
	positions = RandomChoice[Range[0,OptionValue[maxSamples]], nTraces*NReflections];
	magnitudes = RandomReal[{-maxAbsRef, maxAbsRef}, nTraces*NReflections];
	thread = MapThread[SparseArray[(#1+1)->#2, OptionValue[maxSamples]+1]&,{positions, magnitudes}];
	Plus@@@Partition[thread, NReflections]
]


Options[rpWedgeCreateDispersiveTraces]={seed->$RandomSeed,sampleRate -> $SampleRate, maxSamples -> $MaxSamples}
rpWedgeCreateDispersiveTraces[NReflections_, nTraces_, fwav_, maxAbsRef_:0.2, OptionsPattern[]]:=
With[{fr=Subdivide[0., 1/(2.*OptionValue[sampleRate]), OptionValue[maxSamples]]},
	SeedRandom[OptionValue[seed]];
	(*Module[{positions, magnitudes, phasedWlet},
			phasedWlet = Function[{pos, amp}, amp(#/Max@Abs@#)&@rpWedgeInvFourier@rpWedgeSymConjAr@N@Chop@(Exp@(2*I*Pi*fr*pos)rpWedgeWavelet[fr,fwav])];
			positions = RandomChoice[Subdivide[OptionValue[sampleRate]*OptionValue[maxSamples],OptionValue[maxSamples]], nTraces*NReflections];
			magnitudes = RandomReal[{-maxAbsRef, maxAbsRef}, nTraces*NReflections];
			Plus@@@Partition[Threshold[MapThread[phasedWlet,{positions, magnitudes}]//Chop//N, 0.00001][[;;,OptionValue[maxSamples]+1;;]], NReflections]
		]*)
	Module[{positions, magnitudes, Wlet, thread},
			Wlet = (#/Max@Abs@#)&@rpWedgeInvFourier@rpWedgeSymConjAr@N@Chop@rpWedgeWavelet[fr,fwav];
			positions = RandomChoice[Subdivide[OptionValue[maxSamples],OptionValue[maxSamples]], nTraces*NReflections];
			magnitudes = RandomReal[{-maxAbsRef, maxAbsRef}, nTraces*NReflections];
			thread = MapThread[#2 RotateRight[Wlet,#1]&,{positions, magnitudes}]//Chop//N;
			Plus@@@Partition[Threshold[thread, 0.00001][[;;,OptionValue[maxSamples]+1;;]], NReflections]
		]
]


Options[rpWedgeCreateNoise]={seed->$RandomSeed, sampleRate -> $SampleRate, maxSamples -> $MaxSamples}
rpWedgeCreateNoise[nTraces_, fLimit_, maxLevel_:0.01, OptionsPattern[]]:=
	Module[{fourierSeries, bandLimitNoise},
			SeedRandom[OptionValue[seed]];
			fourierSeries = RandomReal[{0.0001, 1}, 7*fLimit];
			bandLimitNoise = PadRight[fourierSeries, OptionValue[maxSamples] + IntegerPart[3*7*fLimit] + 1];
			(maxLevel #/Max@Abs@#)&@(rpWedgeInvFourier@rpWedgeSymConjAr[bandLimitNoise])[[OptionValue[maxSamples] + 2*IntegerPart[3*7*fLimit] + 3 ;;]]//Chop//N
		]


Options[rpThinLayerQ]={sampleRate -> $SampleRate}
rpThinLayerQ[f0_, trace_?VectorQ]:= With[{thinLength = Ceiling[1/f0/OptionValue[sampleRate]]},
	With[{padtrace = PadRight[trace, thinLength]},
		padtrace(*Pick[#, Unitize[#], 1]&/@Partition[padtrace, thinLength,1]*)
	]
]


Options[rpInitiateThinLayerDictionary]={sampleRate->$SampleRate, maxSamples->$MaxSamples, threshold->0.001, scaleFactors->{1., 0.9}, scaleType->"unbiased"};
rpInitiateThinLayerDictionary[maxT_, tSteps_, freq_, sSteps_Integer:1, opts:OptionsPattern[]]:=
With[{w=rpWedgeElastWedgeRef[maxT, tSteps, freq, FilterRules[{opts},Options[rpWedgeElastWedgeRef]]], thres=OptionValue[threshold]},
	With[{rot=Table[RotateRight[i,j][[OptionValue[maxSamples]+1;;]],{j,1,(Length@First@w-1)/2, sSteps},{i,w}]},
		SparseArray@Transpose@Threshold[Flatten[rot,1],thres]
	]
]


(* ::Input::Initialization:: *)
End[]


(* ::Input::Initialization:: *)
EndPackage[]
