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
highlightElementsTable::usage = "";
highlightMatrixElements::usage = "";
calculateMatrixDims::usage = "";
transformToMatrix::usage = "";
getRandomSystem::usage = "";
exerciseSystemToMatrix::usage = "";

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

highlightElementsTable[eqs_] := Module[{rowCount,colCount,incognite,coefficienti,termineNoto},
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
	Style[TableForm[eqsData,TableSpacing->{3,1}, TableAlignments->{Right,Left}],30]
];

highlightMatrixElements[matrix_]:= Module[{lastCol,editedMatrix},
	lastCol = 1;
	For[i=1,i<=Length[matrix],i++,
		If[Length[matrix[[i]]] > lastCol,lastCol = Length[matrix[[i]]]];
	];
	editedMatrix = MapAt[Style[#,RGBColor[0.13,0.52,0.96]]&,matrix,{{Range[1,Length[matrix]],Range[1,lastCol-1]}}];
	editedMatrix = MapAt[Style[#,RGBColor[0.14,0.61,0.14]]&,editedMatrix,{{Range[1,Length[matrix]],lastCol}}];
	Return[editedMatrix];
];

calculateMatrixDims[system_] := Module[{rows,cols,lhsParts},
	rows = Length[system];
	lhsParts = Level[#,1][[1]] & /@ system;
	cols = Length[Variables[lhsParts]]+1;
	Return[{rows,cols}];
];

getRandomSystem[] := Module[{systems},
	systems = {
		{3x-2y==-1, 4x-5y==-2},
		{5y+x==3, 2x-4y==-8},
		{2x-5y==7, x-3y==1},
		{2x-y==0, x+3y==1},
		{z+3x-2y==0, z+x-y==0, 4x+2y-3z==5},
		{x+y-z==6, y+x==3, x+z==0},
		{3x+y+z==3, 6x-2y+z==1, 3y+3x+3z==7},
		{x+3z+2y==1, 3x+4y+6z==3, 5y-3z+10x==-4}
	};
	Return[RandomChoice[systems]];
];

transformToMatrix[eqs_,withTerms_:False] := Module[{system,matrix,terms,row,incognite,rules,i,j,key},
	If[withTerms,
		system = Level[#,1][[1]] & /@ eqs;
		terms = Level[#,1][[2]] & /@ eqs,
		system = eqs
	];
	incognite = Sort[Variables[system]];
	rules = CoefficientRules[system,incognite];
	matrix = {};
	For[i=1,i<=Length[rules],i++,
		row = {};
		For[j=1,j<=Length[incognite],j++,
			key = Normal[SparseArray[{j->1},Length[incognite]]];
			If[MemberQ[Keys[rules[[i]]],key],
				AppendTo[row,Lookup[rules[[i]],Key[key]]],
				AppendTo[row,0]
			];
		];
		If[withTerms,AppendTo[row,terms[[i]]]];
		AppendTo[matrix,row];
	];
	Return[matrix];
];

exerciseSystemToMatrix[] := Module[{system,matrix,dims,rowCount,colCount,inputMatrix,inputMatrixShown,checkButton},
	system = getRandomSystem[];
	matrix = transformToMatrix[system,True];
	dims = calculateMatrixDims[system];
	rowCount = dims[[1]];
	colCount = dims[[2]];
	inputMatrix = ConstantArray[0,rowCount*colCount];
	inputMatrixShown = Partition[Table[With[{i=i},InputField[Dynamic[inputMatrix[[i]]],FieldSize->{1,1}]],{i,rowCount*colCount}],colCount];
	checkButton = Button["Verifica!", If[MatchQ[matrix,ArrayReshape[inputMatrix,dims]],MessageDialog["CORRETTO!"],MessageDialog["SBAGLIATO! Riprova!"]]];
	Grid[{{displayEquationSystem[system],inputMatrixShown//MatrixForm,checkButton}}]
];

End[];
Protect["GaussLinearSystemsPackage`*"]
EndPackage[];
