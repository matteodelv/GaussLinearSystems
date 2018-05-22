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
displayEquationSystem::usage = "Permette la visualizzazione di una lista di equazioni sotto forma di sistema";
highlightElementsTable::usage = "Partendo da un sistema di equazioni, genera una tabella contenente coefficienti, incognite e termini noti, opportunamente evidenziati e colorati";
highlightMatrixElements::usage = "Applica gli stili ad una matrice per evidenziare coefficienti e termini noti";
calculateMatrixDims::usage = "Restituisce la dimensioni (righe, colonne) di una matrice";
transformToMatrix::usage = "Partendo da un sistema di equazioni, restituisce la matrice associata";
getRandomSystem::usage = "Restituisce un sistema di equazioni random tra quelli presenti nel dataset";
exerciseSystemToMatrix::usage = "Mostra l'esercizio in cui viene visualizzato un sistema di equazioni random e l'utente deve trasformarlo nella matrice associata";
exerciseReducedMatrixToSystem::usage = "Mostra l'esercizio in un viene visualizzato un sistema di equazioni ridotto con Gauss e l'utente deve ricostruire il sistema ridotto associato";
exerciseTriangularizeMatrix::usage = "Mostra l'esercizio in cui viene visualizzata una matrice (relativa ad un sistema random) e l'utente deve applicare il Metodo di Gauss per inserire la matrice ridotta associata";
exerciseFinalGauss::usage = "Mostra l'esercizio finale in cui l'utente inserisce un sistema di equazioni e applica il Metodo di Gauss passo per passo";
fattorizzazioneLU::usage = "Implementazione della Fattorizzazione LU con pivoting parziale. Partendo da una matrice quadrata A ed un vettore b, ritorna le matrici L ed U, la matrice di permutazione P ed il vettore b' permutato";
diagonalMatrixQuestion::usage="Domanda: Trovare la diagonale di una matrice";
rowColumnQuestion::usage="Domanda: Trovare quante righe/colonne possiede la matrice indicata";
indexQuestion::usage="Domanda: Indentificare un elemento di una matrice, dato l'indice o il valore. ";
identityMatrixQuestion::usage="Domanda: Viene richiesto dove la matrice identit\[AGrave] presenta degli elementi diversi da zero(Diagonale)";
haveDiagonalQuestion::usage="Domanda: Viene richiesto quale tra le matrice fornite possiede una diagonale ";
questionsExercise::usage="Gestisce la casualit\[AGrave] delle domande, la correttezza della risposta fornita e la relativa stampa sul Notebook";
checkButtonFunction::usage = "";

oneElementList::usage = "";
stringInputToSystem::usage = "";

t::usage = "";
x::usage = "";
y::usage = "";
z::usage = "";

Begin["`Private`"];

Off[Solve::svars];
Off[General::shdw];
SetDirectory[NotebookDirectory[]];

(*Funzione che permette di mostrare graficamente la soluzione 
di un sistema di due equazioni di primo grado in due incognite *)
plotLinearSystem2[eq1_, eq2_] := Module[{solutions, points},
	solutions = Solve[{eq1,eq2}];
	points = {Blue, PointSize[Large], Point[{x,y}/.solutions]};
	plotEq1 = Plot[y /. Solve[eq1], {x,-10,10},AspectRatio->1,PlotStyle->Red,PlotLegends->{eq1//TraditionalForm}];
	plotEq2 = Plot[y /. Solve[eq2], {x,-10,10},AspectRatio->1,PlotStyle->Orange,PlotLegends->{eq2//TraditionalForm}];
	Show[plotEq1, plotEq2, Graphics[{points}], ImageSize->400,Background->White]
];

(* Funzione che permette di mostrare graficamente la soluzione 
di un sistema di tre equazioni di primo grado in tre incognite *)
plotLinearSystem3[eq1_,eq2_,eq3_] := Module[{sols, points},
	sols = Solve[{eq1,eq2,eq3}];
	points = {Red, PointSize[Large], Point[{x,y,z} /. sols]};
	plotEq1 = Plot3D[z /. Solve[eq1], {x,-10,10}, {y,-10,10},AspectRatio->1,Mesh->None,PlotStyle->RGBColor[1.0,1.0,0,0.4],PlotLegends->{eq1//TraditionalForm}];
	plotEq2 = Plot3D[z /. Solve[eq2], {x,-10,10}, {y,-10,10},AspectRatio->1,Mesh->None,PlotStyle->RGBColor[0.6,0.2,1.0,0.4],PlotLegends->{eq2//TraditionalForm}];
	plotEq3 = Plot3D[z /. Solve[eq3], {x,-10,10}, {y,-10,10},AspectRatio->1,Mesh->None,PlotStyle->RGBColor[0.99,0.47,1.0,0.4],PlotLegends->{eq3//TraditionalForm}];
	Show[plotEq1,plotEq2,plotEq3,Graphics3D[{points}],ImageSize->400]
];

(*Funzione che data in input una lista di equazione la visualizza sotto forma di sistema*)
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

(* Funzione che permette di evidenziare gli elementi di una matrice *)
highlightMatrixElements[matrix_]:= Module[{lastCol,editedMatrix},
	lastCol = 1;
	For[i=1,i<=Length[matrix],i++,
		If[Length[matrix[[i]]] > lastCol,lastCol = Length[matrix[[i]]]];
	];
	editedMatrix = MapAt[Style[#,RGBColor[0.13,0.52,0.96]]&,matrix,{{Range[1,Length[matrix]],Range[1,lastCol-1]}}];
	editedMatrix = MapAt[Style[#,RGBColor[0.14,0.61,0.14]]&,editedMatrix,{{Range[1,Length[matrix]],lastCol}}];
	Return[editedMatrix];
];

(* Funzione che restituisce la dimensioni (righe, colonne) di una matrice *)
calculateMatrixDims[system_] := Module[{rows,cols,lhsParts},
	rows = Length[system];
	lhsParts = Level[#,1][[1]] & /@ system;
	cols = Length[Variables[lhsParts]]+1;
	Return[{rows,cols}];
];

(* Restituisce un sistema random tra quelli presenti. Funzione utilizzata negli esercizi*)
getRandomSystem[solvable_:True] := Module[{systems},
	systems = {
		{3x-2y==-1, 4x-5y==-2} (* 2 equazioni in 2 incognite *),
		{5y+x==3, 2x-4y==-8},
		{2x-5y==7, x-3y==1},
		{2x-y==0, x+3y==1},
		{5x+y==20, 5x+7y==20},
		{2x-3y==4, -2x+4y==1},
		{2x-y==7, 4x+3y==4},
		{1/3x+4y==5, -x+1/2y==-5/2},
		{x-1/2y==5/3, 3/2x-3/8y==1},
		{z+3x-2y==0, z+x-y==0, 4x+2y-3z==5} (* 3 equazioni in 3 incognite *),
		{x+y-z==6, y+x==3, x+z==0},
		{3x+y+z==3, 6x-2y+z==1, 3y+3x+3z==7},
		{x+3z+2y==1, 3x+4y+6z==3, 5y-3z+10x==-4},
		{2x+3y-z==0, x-y+z==1, 3x+2y+4z==-3},
		{x+z==1, y+1/2z==-1, 2x+z==0},
		{-x+y-z==1, 10x-5y+10z==-3, 2x+y+z==2},
		{3x-y+2z==10, 6x+4z-y==17, x-2z+2y==-5},
		{y+z-t==1, x-2y+t==1, 3x+2y-z-t==0, x-z==-2} (* 4 equazioni in 4 incognite *),
		{2t+x+3y+5z==4, 4t+2y+8z==8, t+2x+2y+3z==5, x+y+z==4},
		{2t-3x+4z==-1, x+2y-t==2, 3y+2z+t==4, x+y+t==0}
	};
	If[Not[solvable],
		systems = Join[systems, {
			{1/2x-y+2t==0, y+2z-2t==1, x+y+6z-2t==3, y-z+3t==0} (* impossibile *),
			{3x+2y-z==0, x+y+z==3} (* indeterminato *),
			{x-2y+z==1, -t+3y+z==2, 7t+7z==7, y+2z==3} (* impossibile *),
			{t-x+y==1, 2t-y==0, -t+z==-1} (* indeterminato *),
			{2x+3y+t==6, -9y+2z-t==3, -z+7t==5} (* indeterminato *),
			{x+y==0, x+y==1} (* impossibile *),
			{x-y+z==2, -x+y+z==1} (* indeterminato *)
		}];
	];
	Return[RandomChoice[systems]];
];
(*Funzione che prende in input un sistema di equazioni e restituisce la matrice associata*)
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
(*Funzione esercizio: che dato un sistema di equazioni random l'utente deve trasformarlo nella matrice associata*)
exerciseSystemToMatrix[inputSystem_:Null,result_:0] := DynamicModule[
	{dims, rowCount, colCount, checkButtonSystemToMatrix, restartButton, gridOptions, system, inputMatrixShown, matrix, inputMatrix,
	okColor=RGBColor[0,1,0,0.4], wrongColor=RGBColor[1,0,0,0.4], shown=False, dialogImage, dialogText, output},
	(*viene scelto un sistema random tra quelli presenti nel sistema*)
	If[SameQ[inputSystem,Null], system = getRandomSystem[False], system = inputSystem];
	(*calcola la dimensione della matrice*)
	{rowCount,colCount} = calculateMatrixDims[system];
	matrix = transformToMatrix[system,True];
	matrix = Flatten[matrix];
	inputMatrix = ConstantArray[0,rowCount*colCount];
	inputMatrixShown = Partition[
		Table[With[{i=i},
			Dynamic[Framed[
				InputField[Dynamic[inputMatrix[[i]]], FieldSize->{2,1}, Appearance->"Frameless",
					Background->If[shown,
						If[SameQ[inputMatrix[[i]],matrix[[i]]], okColor, wrongColor],
						White]
					],
					FrameStyle->If[shown,
						If[SameQ[inputMatrix[[i]],matrix[[i]]], okColor, wrongColor],
						Automatic],
					FrameMargins->None
				]]
			], {i,rowCount*colCount}
		], colCount];
	(*Button check: controlla che la matrice inserita corrisponda alla soluzione*)
	checkButtonSystemToMatrix = Button[Style["Verifica!",24],
		Clear[dialogImage, dialogText];
		Print["checkButtonSystemToMatrix"];
		shown=False;
		If[MatchQ[matrix,inputMatrix],
			Print["IF"];
			If[Not[SameQ[inputMatrix, Null]], result = 2; Print[result]];
			dialogImage = Import["images/checkmark.png"];
			dialogText = Style["CORRETTO. Bravo!",20,RGBColor[0.14,0.61,0.14]],
			If[Not[SameQ[inputMatrix, Null]], result = 1; Print[result]];
			dialogImage = Import["images/error.png"];
			dialogText = Style["SBAGLIATO. Riprova!",20,Red];
			shown=True;
		];
		MessageDialog[Column[{dialogImage,dialogText}, Spacings->{2,4}, Alignment->Center,
			Frame->All, FrameStyle->RGBColor[0,0,0,0], ItemSize->Fit]],
		ImageSize->{200,50}];
	(*Button ricomincia. Viene rivalutato il notebook *)
	If[SameQ[inputSystem, Null],
		restartButton = Button[Style["Ricomincia!",24],
			shown=False;
			NotebookFind[EvaluationNotebook[], "exerciseSystemToMatrixCellTag",All,CellTags];
			SelectionEvaluate[EvaluationNotebook[]],
			ImageSize->{200,50}];
	];
	gridOptions = {
		Frame->All,
		FrameStyle->RGBColor[0.94,0.94,0.94],
		ItemStyle->Directive[FontFamily->"Roboto Condensed", FontSize->28],
		Spacings->{10,3},
		ItemSize->Fit,
		Alignment->{{Right,Left}}
	};
	output = {
		{displayEquationSystem[system],inputMatrixShown//MatrixForm},
		{checkButtonSystemToMatrix,If[SameQ[inputSystem,Null], restartButton,SpanFromLeft]}
	};
	If[Not[SameQ[inputSystem,Null]],
		Grid[output, gridOptions],
		Grid[output, gridOptions]
	]
];
SetAttributes[exerciseSystemToMatrix,HoldRest];

exerciseReducedMatrixToSystem[system_:Null,result_:0] := DynamicModule[
	{reducedMatrix, inputSystem, inputSystemShown, randomSystem, A, b, U, L, P, newB, matrix, rows, checkButtonReduced, cols, prova, output,
	okColor=RGBColor[0,1,0,0.4], wrongColor = RGBColor[1,0,0,0.4], shown = False, backup,provaMatrice, restartButton, gridOptions, dialogImage, dialogText},
	
	If[SameQ[system, Null], randomSystem = getRandomSystem[], randomSystem = system];
	{rows, cols} = calculateMatrixDims[randomSystem];
	matrix = transformToMatrix[randomSystem,True];
	
	A = Drop[matrix, None, -1];
	b = Take[matrix, All, -1];
	{L,U,P,newB} = fattorizzazioneLU[A,b];
	reducedMatrix = Join[U, ArrayReshape[newB,{rows,1}],2];
	prova = ConstantArray[Null, rows];
	provaMatrice = ArrayReshape[ConstantArray[0,rows*cols],{rows, cols}];
	
	inputSystem = ConstantArray[Null, rows];
	inputSystemShown = Partition[
			Table[With[{i=i},
				Dynamic[Framed[
					InputField[Dynamic[inputSystem[[i]]], String,FieldSize->{15,1},Appearance->"Frameless",
						Background->If[shown,
							If[SameQ[reducedMatrix[[i]],provaMatrice[[i]]],RGBColor[0,1,0,0.4],RGBColor[1,0,0,0.4]],
							White]
						],
						FrameStyle->If[shown,
							If[SameQ[reducedMatrix[[i]],provaMatrice[[i]]],RGBColor[0,1,0,0.4],RGBColor[1,0,0,0.4]],
							Automatic],
						FrameMargins->None
					]]
				], {i,rows}
			],1];
	checkButtonReduced = Button[Style["Verifica!",24],
		If[Not[MemberQ[inputSystem, Null]],
			Clear[dialogImage, dialogText];
			prova = StringReplace[#, "=" -> "\[Equal]"] & @ inputSystem;
			prova = ToExpression[#] & @ prova;
			provaMatrice = transformToMatrix[prova, True];
			If[MatchQ[reducedMatrix,provaMatrice],
				If[Not[SameQ[system, Null]], result = 4];
				dialogImage = Import["images/checkmark.png"];
				dialogText = Style["CORRETTO. Bravo!",20,RGBColor[0.14,0.61,0.14]],
				If[Not[SameQ[system, Null]], result = 3];
				dialogImage = Import["images/error.png"];
				dialogText = Style["SBAGLIATO. Riprova!",20,Red];
			],
			If[Not[SameQ[system, Null]], result = 3];
			dialogImage = Import["images/error.png"];
			dialogText = Style["SBAGLIATO. Riprova!",20,Red];
		];
		MessageDialog[Column[{dialogImage,dialogText}, Spacings->{2,4}, Alignment->Center,
				Frame->All, FrameStyle->RGBColor[0,0,0,0], ItemSize->Fit]];
		shown = True,
		ImageSize->{200,50}
	];
	If[SameQ[system, Null],
		restartButton = Button[Style["Ricomincia!",24],
			shown=False;
			NotebookFind[EvaluationNotebook[], "exerciseReducedMatrixToSystemCellTag",All,CellTags];
			SelectionEvaluate[EvaluationNotebook[]],
			ImageSize->{200,50}];	
	];
	gridOptions = {
		Frame->All,
		FrameStyle->RGBColor[0.94,0.94,0.94],
		ItemStyle->Directive[FontFamily->"Roboto Condensed", FontSize->28],
		Spacings->{10,3},
		ItemSize->Fit,
		Alignment->{{Right,Left}}
	};
	output = {
		{highlightMatrixElements[reducedMatrix]//MatrixForm,inputSystemShown//MatrixForm},
		{checkButtonReduced,If[SameQ[system, Null],restartButton,SpanFromLeft]}
	};
	If[Not[SameQ[system, Null]],
		Grid[output, gridOptions],
		Grid[output, gridOptions]
	]
];
SetAttributes[exerciseReducedMatrixToSystem, HoldRest];

exerciseTriangularizeMatrix[system_, finalExercise_:False, result_:0] := DynamicModule[
	{systemMatrix, rows, cols, A, b, L, U, P, newB, UFlattened, shown=False, okColor=RGBColor[0,1,0,0.4], 
	wrongColor=RGBColor[1,0,0,0.4], inputMatrix, inputMatrixShown, dialogImage, dialogText, checkButtonTriangular, restartButton, gridOptions, output},
	
	systemMatrix = transformToMatrix[system,True];
	{rows,cols} = calculateMatrixDims[system];
	A = Drop[systemMatrix, None, -1];
	b = Take[systemMatrix, All, -1];
	{L,U,P,newB} = fattorizzazioneLU[A,b];
	
	If[Not[MatchQ[{L,U,P,newB}, ConstantArray[Null,4]]]
		(* Gestire caso sistema indeterminato o impossibile *),
		UFlattened = Flatten[Join[U,ArrayReshape[newB,{rows,1}],2]];
		inputMatrix = ConstantArray[0, rows*cols];
		inputMatrixShown = Partition[
			Table[With[{i=i},
				Dynamic[Framed[
					InputField[Dynamic[inputMatrix[[i]]],FieldSize->{2,1},Appearance->"Frameless",
						Background->If[shown,
							If[SameQ[inputMatrix[[i]],UFlattened[[i]]],okColor,wrongColor],
							White]
						],
						FrameStyle->If[shown,
							If[SameQ[inputMatrix[[i]],UFlattened[[i]]],okColor,wrongColor],
							Automatic],
						FrameMargins->None
					]]
				], {i,rows*cols}
			],cols];
		For[i=2,i<=Length[inputMatrixShown],i++,
			inputMatrixShown[[i]][[1;;i-1]]=0;
		];
		checkButtonTriangular = Button[Style["Verifica!",24],
			Clear[dialogImage, dialogText];
			shown=False;
			If[MatchQ[UFlattened,inputMatrix],
				If[finalExercise, result = 3];
				dialogImage = Import["images/checkmark.png"];
				dialogText = Style["CORRETTO. Bravo!",20,RGBColor[0.14,0.61,0.14]],
				If[finalExercise, result = 2];
				dialogImage = Import["images/error.png"];
				dialogText = Style["SBAGLIATO. Riprova!",20,Red];
				shown=True;
			];
			MessageDialog[Column[{dialogImage,dialogText}, Spacings->{2,4}, Alignment->Center,
				Frame->All, FrameStyle->RGBColor[0,0,0,0], ItemSize->Fit]],
			ImageSize->{200,50}];
		If[Not[finalExercise],
			restartButton = Button[Style["Ricomincia!",24],
				shown=False;
				NotebookFind[EvaluationNotebook[], "exerciseTriangolarizzazioneCellTag",All,CellTags];
				SelectionEvaluate[EvaluationNotebook[]],
				ImageSize->{200,50}];
		];
		gridOptions = {
			Frame->All,
			FrameStyle->RGBColor[0.94,0.94,0.94],
			ItemStyle->Directive[FontFamily->"Roboto Condensed", FontSize->28],
			Spacings->{10,3},
			ItemSize->Fit,
			Alignment->{{Right,Left}}
		};
		
		output = {
			{highlightMatrixElements[systemMatrix]//MatrixForm,inputMatrixShown//MatrixForm},
			{checkButtonTriangular,If[Not[finalExercise],restartButton,SpanFromLeft]}
		};
		If[Not[finalExercise],
			Grid[output,gridOptions],
			Grid[output, gridOptions]
		]
	]
];
SetAttributes[exerciseTriangularizeMatrix, HoldRest];

(*Fattorizzazione LU: *)
fattorizzazioneLU[A_,b_] := Module[{L,U,P,bPerm,matriceEdited,n,pivot,candidatePivot,subColonna,pivotIndex,lambda,i,k},
	If[SameQ[Det[A],0]||Not[Equal[Length[A],Length[A[[All,1]]]]],Return[{Null,Null,Null,Null}]];

	n = Length[A];
	matriceEdited = A;
	P = IdentityMatrix[n];
	L = ArrayReshape[ConstantArray[0,n*n],{n, n}];
	For[i=1,i<n,i++,
		subColonna = matriceEdited[[All,i]][[i;;]];
		candidatePivot = Max[Abs[subColonna]];
		pivotIndex = FirstPosition[Abs[matriceEdited[[All,i]]][[i;;]],candidatePivot][[1]];
		pivotIndex = pivotIndex+i-1;
		pivot = matriceEdited[[pivotIndex,i]];
		If[Not[Equal[i,pivotIndex]],
			matriceEdited = Permute[matriceEdited,Cycles[{{i,pivotIndex}}]];
			P = Permute[P, Cycles[{Flatten[{i,pivotIndex}]}]];
		];
		For[k=i+1,k<=n,k++,
			If[Not[Equal[matriceEdited[[k,i]],0]],
				lambda = (matriceEdited[[k,i]]/pivot);
				L[[k,i]] = lambda;
				matriceEdited[[k]] = matriceEdited[[k]]-lambda*matriceEdited[[i]];
			];
		];
		If[Not[Equal[i,pivotIndex]]&&i>=2,
			L[[i;;n,1;;i-1]]=Permute[L[[i;;n,1;;i-1]],Cycles[{{1,pivotIndex-i+1}}]];
		];
	];
	U = matriceEdited;
	L=L+IdentityMatrix[n];
	bPerm = Inverse[L].P.b;
	Return[{L,U,P,bPerm}];
];

(* Data una matrice generata random viene richiesto all'utente di selezionarne la diagonale*)
diagonalMatrixQuestion[]:= Module[{index, randomElements,matrix, answers,solution,wrongAnswers,text, gridOptions},
	text = "Quale tra queste \[EGrave] la diagonale della seguente matrice ?";
	(*Generazione grandezza matrice, elementi, risposte errate e soluzioni*)
	index = RandomInteger[{3,5}];
	randomElements = RandomSample[Range[-20,20],(index*index)];
	matrix = ArrayReshape[randomElements,{index,index}];
	(*La soluzione corrisponde alla diagonale della matrice*)
	solution = Diagonal[matrix];
	wrongAnswers = RandomSample[Tuples[DeleteCases[randomElements,solution],index],3];
	(*Unione della risposta esatta con le risposte sbagliate*)
	answers = Join[{solution},wrongAnswers];
	Return[{text,MatrixForm/@answers,solution//MatrixForm,matrix}]			
];
(*Funzione che data una matrice di grandezza casuale richiede all'utente quale delle due indica la colonna o la riga*)
rowColumnQuestion[]:= Module[{index,row,column,matrix,answers,text,solution},
	(*Generazione di una matrice non quadrata e della soluzione casuale tra riga e colonna*)
	index = RandomSample[Range[4,6],2];
	row = index[[1]];
	column = index[[2]];
	matrix = (row*column);
	answers = {row, column, matrix, (row+column)};
	solution = RandomChoice[{row,column}];
	text = StringReplace["Quante righe ci sono in una matrice NRighe X NColonne ?", 
			{"NRighe"-> ToString[row],"NColonne"-> ToString[column]}];
	If[MatchQ[solution,column],
		text = StringReplace[text, "righe"->"colonne"]
	];
	Return[{text,answers,solution, Null}]
];

(*Domanda: Indentificare un elemento di una matrice, dato l'indice o il valore*)
indexQuestion[]:=Module[{index,row,column,position,matrixRandom,matrix,solution,text,element, 
						answers,wrongAnswersRow,wrongAnswersColumn,wrongAnswers, ind, subscriptIndex},
	(*Generazione matrice di grandezza casuale e scelta casuale di una posizione*)
	index = RandomSample[Range[4,6],2];
	row = RandomInteger[{1,index[[1]]}];
	column = RandomInteger[{1,index[[2]]}];
	ind = StringReplace["xy",{"x"->ToString[row],"y"->ToString[column]}];
	position = Subscript["a","index"]/.{"index"->ind};
	(*Viene popolata casualmente la matrice*)
	matrixRandom = RandomSample[Range[-15,15],(index[[1]]*index[[2]])];
	matrix = ArrayReshape[matrixRandom,{index[[1]],index[[2]]}];
	(*Casualit\[AGrave] della domanda tra elemento teorico o pratico*)
	element = RandomChoice[{position//TraditionalForm,matrix[[row,column]]}];
	text= StringReplace["Nella matrice l'elemento x a cosa corrisponde ?","x"-> ToString[element]];
	If[MatchQ[element,position//TraditionalForm],
		(*Generazione random delle soluzioni errate*)
		solution = matrix[[row,column]];
		wrongAnswers = RandomSample[Tuples[DeleteCases[matrixRandom,solution],1],3],
			solution = position;
			wrongAnswersRow = RandomSample[Tuples[Range[1,index[[1]]],1],3];
			wrongAnswersColumn = RandomSample[Tuples[DeleteCases[Range[1,index[[2]]],column],1],3];
			wrongAnswers = Join[wrongAnswersRow,wrongAnswersColumn,2];
			ind = StringReplace["xy",{"x"->ToString[#[[1]]],"y"->ToString[#[[2]]]}] & /@ wrongAnswers;
			wrongAnswers = Subscript["a","index"]/.{"index"->#} & /@ ind	
	];
	answers = Flatten[{solution, wrongAnswers}];
	Return[{text,answers,solution,matrix}]
];

(*Domanda che richiede all'utente la caratteristia chiave delle matrici identit\[AGrave]*)
identityMatrixQuestion[]:=Module[{text,solution,wrongAnswers,answers},
	text = "La matrice identit\[AGrave] presenta degli elementi diversi da zero:";
	solution = "Nella diagonale";
	wrongAnswers = {"Nella prima riga","Nella prima colonna","Indifferente, purch\[EGrave] siano 1"};
	(*Vengono unite la soluzione e le risposte errate*)
	answers = Flatten[{solution,wrongAnswers}];
	Return[{text,answers,solution,Null}]
];

(*Domanda che richiede all'utente quale matrice possiede una diagonale *)
haveDiagonalQuestion[]:= Module[{text,solution,wrongAnswers,answers},
	text = "Quale delle seguenti matrici possiede una DIAGONALE ?";
	solution = "Matrice quadrata";
	wrongAnswers = {"Matrice rettangolare","Entrambe","Nessuna"};
	(*Vengono unite la soluzione e le risposte errate*)
	answers = Flatten[{solution,wrongAnswers}];
	Return[{text,answers,solution,Null}]
];

(*Funzione che sceglie random le domande a risposta multipla e ne gestisce la relativa stampa sul Notebook*)
questionsExercise[]:=DynamicModule[{question,text,answers,solution,matrix,radioAnswers,choice,checkButton,
									restartButton,gridOptions, appearance,dialogImage,dialogText},
	(*Viene scelta la domanda random*)
	question = RandomChoice[{
		diagonalMatrixQuestion,
		haveDiagonalQuestion,
		rowColumnQuestion,
		identityMatrixQuestion,
		indexQuestion
	}];
	{text,answers,solution,matrix} = question[];
	answers = RandomSample[answers];
	(*RadioButtonBar con le risposte e relativo appeareance in base alla domanda*)
	If[MatchQ[question,diagonalMatrixQuestion], 
		appearance="Horizontal",
		appearance="Vertical"
	];
	radioAnswers= Dynamic[RadioButtonBar[Dynamic[choice],answers, Appearance->appearance]];
	(*Controllo della scelta effettuata*)
	checkButton = Button[Style["Verifica!",24],
		If[MatchQ[choice,solution],	
			dialogImage = Import["images/checkmark.png"];
			dialogText = Style["CORRETTO. Bravo!",20,RGBColor[0.14,0.61,0.14]],
			dialogImage = Import["images/error.png"];
			dialogText = Style["SBAGLIATO. Riprova!",20,Red];
		];
		MessageDialog[Column[{dialogImage,dialogText}, Spacings->{2,4}, Alignment->Center,
			Frame->All, FrameStyle->RGBColor[0,0,0,0], ItemSize->Fit]],
		ImageSize->{200,50}];
	(*Generazione nuova domanda con rivalutazione del notebook nel tag specificato*)
	restartButton = Button[Style["Nuova domanda",24],
		NotebookFind[EvaluationNotebook[], "questionExerciseTag",All,CellTags];
		SelectionEvaluate[EvaluationNotebook[]],
		ImageSize->{300,50}];
	(*Stampa della domanda con personalizzazione della grid di Output*)
	gridOptions = {
		Frame->All,
		FrameStyle->RGBColor[0.94,0.94,0.94],
		ItemStyle->Directive[FontFamily->"Roboto Condensed", FontSize->28],
		Spacings->{10,3},
		ItemSize->Fit,
		Alignment->{{Left,Left}}
	};
	If[SameQ[matrix, Null],
		grid = Grid[{{text,SpanFromLeft},{radioAnswers,SpanFromLeft},{restartButton,checkButton}}, gridOptions],
		grid = Grid[{{text,SpanFromLeft},{radioAnswers,matrix//MatrixForm},{restartButton,checkButton}}, gridOptions]
	];
	grid
];

oneElementList[list_, element_] := Module[{},
	If[Not[SameQ[list, {}]], list = {}];
	list = Insert[list, element, 1];
];
SetAttributes[oneElementList, HoldAll];

stringInputToSystem[list_] := Module[{},
	list = StringReplace[#,"="->"\[Equal]"]& @ list;
	list = ToExpression @ Flatten @ StringSplit[#, ","] & @ list;
];
SetAttributes[stringInputToSystem, HoldAll];

exerciseFinalGauss[] := DynamicModule[
	{inputList = {}, startButton, restartButton, step1Shown = False, output, step2Shown = False, step3Shown = False, step4Shown = False, step=0},
	
	startButton = Button["Inizia!", stringInputToSystem[inputList];step=1(*;step1Shown=True*)];
	restartButton = Button["Ricomincia", inputList = {};step=0(*;step1Shown=False*)];
	(*
	If[step\[Equal]0,
		Print["Step = 0"];
		output = Grid[{
			{InputField[Dynamic[Null, oneElementList[inputList,#]&], String, FieldSize\[Rule]{40,2}],SpanFromLeft},
			{startButton,restartButton}
		}];
		Print[output];
	];
	If[step\[Equal]1,
		Print["Step = 1"];
		Print["input = ", input];
		output = Grid[{
			{InputField[Dynamic[Null, oneElementList[inputList,#]&], String, FieldSize\[Rule]{40,2}],SpanFromLeft},
			{startButton,restartButton},
			{"Passo 1: Trasforma il sistema nella matrice associata", SpanFromLeft},
			{exerciseSystemToMatrix[input, step2Shown],SpanFromLeft}
		}];
		Dynamic[If[step2Shown,Print["CIAO"];exerciseFinalGauss[2, input]]];
		Print[output];
	];
	If[step\[Equal]2,
		Print["Step = 2"];
		Print["input = ", input];
		output = Grid[{
			{InputField[Dynamic[Null, oneElementList[inputList,#]&], String, FieldSize\[Rule]{40,2}],SpanFromLeft},
			{startButton,restartButton},
			{"Passo 2: Applica il Metodo di Gauss per ridurre la matrice a gradini", SpanFromLeft},
			{exerciseTriangularizeMatrix[input, True, step3Shown],SpanFromLeft}
		}];
		Dynamic[If[step3Shown, exerciseFinalGauss[3, input]]];
		Print[output];
	];*)
	(*Print[output];*)
	(*
	Dynamic[
		Print["DYNAMIC"];
		If[step4Shown, prova = 4,
			Print["PRIMO ELSE"];
			If[step3Shown, prova = 3,
				Print["SECONDO ELSE"];
				If[step2Shown, prova = 2, Print["TERZO ELSE"]];
			];
		];
	];*)
	
	Grid[{
		{InputField[Dynamic[Null, oneElementList[inputList,#]&], String, FieldSize->{40,2}],SpanFromLeft},
		{startButton,restartButton},
		{Dynamic[
			Switch[step,
			0, "",
			1, "Passo 1: Trasforma il sistema nella matrice associata",
			2, "Passo 2: Applica il Metodo di Gauss per ridurre la matrice a gradini",
			3, "Passo 3: Ricava il sistema associato alla matrice ridotta",
			4, "Passo 4: Scegli la soluzione corretta"]
		],SpanFromLeft},
		{Dynamic[
			Switch[step,
			0, "",
			1, exerciseSystemToMatrix[inputList, step],
			2, exerciseTriangularizeMatrix[inputList, True, step],
			3, exerciseReducedMatrixToSystem[inputList, step],
			4, "SEI ARRIVATO ALL'ULTIMO PASSO!"]
		],SpanFromLeft}
		(*{Dynamic[
			If[step1Shown, exerciseSystemToMatrix[inputList, step2Shown];step1Shown = Not[step2Shown],
				If[step2Shown, exerciseTriangularizeMatrix[inputList, True, step3Shown];step2Shown = Not[step3Shown],
					If[step3Shown, exerciseReducedMatrixToSystem[inputList, step4Shown]; step3Shown = False, ""]
				]
			]
		],SpanFromLeft}*)(*,
		{Dynamic[If[step1Shown, exerciseSystemToMatrix[inputList,step2Shown],""]],SpanFromLeft},
		{Dynamic[If[step2Shown, "Applica il Metodo di Gauss per ridurre la matrice a gradini",""]],SpanFromLeft},
		{Dynamic[If[step2Shown, exerciseTriangularizeMatrix[inputList, True, step3Shown],""]],SpanFromLeft},
		{Dynamic[If[step3Shown, "Applica il Metodo di Gauss per ridurre la matrice a gradini",""]],SpanFromLeft},
		{Dynamic[If[step3Shown, exerciseReducedMatrixToSystem[inputList, step4Shown],""]],SpanFromLeft}*)
	},Frame->All] (* la grid va creata tutta ma le righe devono stare all'interno di dynamic con degli if per farle mostrare a pezzi *)
];

End[];
Protect["GaussLinearSystemsPackage`*"]
EndPackage[];
