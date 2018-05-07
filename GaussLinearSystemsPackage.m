(* ::Package:: *)

(*
* GaussLinearSystemsPackage.m
* Progetto d'esame di Matematica Computazionale (LM Informatica) e Calcolo Numerico e Software Didattico (LM Matematica)
* Anno Accademico 2017/18
*
* Autori: Sara Annese, Laura Cappelli, Matteo Del Vecchio, Matteo Marchesini, Daniele Rossi, Silvia Scibilia
*
* Versione Mathematica utilizzata: 11.2
*)

BeginPackage["GaussLinearSystemsPackage`"];
Unprotect["GaussLinearSystemsPackage`*"];
ClearAll["GaussLinearSystemsPackage`*"];

plotLinearSystem2::usage = "Permette di mostrare graficamente la soluzione di un sistema di due equazioni di primo grado in due incognite";
plotLinearSystem3::usage = "Permette di mostrare graficamente la soluzione di un sistema di tre equazioni di primo grado in tre incognite";
displayEquationSystem::usage = "";
hightlightElementsTable::usage = "";

x::usage = "";
y::usage = "";
z::usage = "";

Begin["`Private`"];

Off[Solve::svars];
Off[General::shdw];
SetDirectory[NotebookDirectory[]];

plotLinearSystem2[eq1_, eq2_] := Module[{solutions, points},
	solutions = Solve[{eq1,eq2},{x,y},Reals];
	points = {Blue, PointSize[Large], Point[{x,y} /. solutions]};
	plotEq1 = Plot[y /. Solve[eq1], {x,-10,10},AspectRatio->1,PlotStyle->Red,PlotLegends->{eq1//TraditionalForm}];
	plotEq2 = Plot[y /. Solve[eq2], {x,-10,10},AspectRatio->1,PlotStyle->Orange,PlotLegends->{eq2//TraditionalForm}];
	Show[plotEq1, plotEq2, Graphics[{points}], ImageSize->400,Background->White]
];

plotLinearSystem3[eq1_,eq2_,eq3_] := Module[{sols, points},
	sols = Solve[{eq1,eq2,eq3}];
	points = {Red, PointSize[Large], Point[{x,y,z} /. sols]};
	plotEq1 = Plot3D[z /. Solve[eq1], {x,-10,10}, {y,-10,10},AspectRatio->1,Mesh->None,PlotStyle->RGBColor[1.0,1.0,0,0.4],PlotLegends->{eq1//TraditionalForm}];
	plotEq2 = Plot3D[z /. Solve[eq2], {x,-10,10}, {y,-10,10},AspectRatio->1,Mesh->None,PlotStyle->RGBColor[0.6,0.2,1.0,0.4],PlotLegends->{eq2//TraditionalForm}];
	plotEq3 = Plot3D[z /. Solve[eq3], {x,-10,10}, {y,-10,10},AspectRatio->1,Mesh->None,PlotStyle->RGBColor[0.99,0.47,1.0,0.4],PlotLegends->{eq3//TraditionalForm}];
	Show[plotEq1,plotEq2,plotEq3,Graphics3D[{points}],ImageSize->400]
];

displayEquationSystem[eqs_] := Module[{eqsFF},
	eqsFF = TraditionalForm /@ HoldForm /@ eqs;
	DisplayForm@RowBox[{StyleBox["{", SpanMaxSize->Infinity], Column[eqsFF, Alignment->Left]}]
];

hightlightElementsTable[eqs_] := Module[{rowCount,colCount,incognite,coefficienti,termineNoto},
	rowCount = Length[eqs];
	colCount=0;
	eqsData = {};
	For[i=1, i<=rowCount,i++,
		eqData = {};
		lhs = Level[eqs[[i]],1][[1]];
		termineNoto = Level[eqs[[i]],1][[2]];
		incognite = Sort[Variables[lhs]];
		coefficienti = Values[CoefficientRules[lhs,incognite]];
		If[Length[incognite]*2 > colCount,colCount=Length[incognite]*2];
		For[k=1,k<=Length[incognite],k++,
			AppendTo[eqData, Style[coefficienti[[k]], RGBColor[0.13,0.52,0.96]]];
			AppendTo[eqData, Style[incognite[[k]], Red]];
		];
		AppendTo[eqData, Style[termineNoto, RGBColor[0.14,0.61,0.14]]];
		eqsData = Append[eqsData, eqData];
	];
	Style[TableForm[eqsData,TableSpacing->{3,1}, TableAlignments->{Right,Left}],Large]
];

End[];
Protect["GaussLinearSystemsPackage`*"]
EndPackage[];
