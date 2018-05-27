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

(* Funzioni principali *)
fattorizzazioneLU::usage = "Implementazione della Fattorizzazione LU con pivoting parziale. Partendo da una matrice quadrata A ed un vettore b, ritorna le matrici L ed U, la matrice di permutazione P ed il vettore b' permutato";
plotLinearSystem2::usage = "Permette di mostrare graficamente la soluzione di un sistema di due equazioni di primo grado in due incognite";
plotLinearSystem3::usage = "Permette di mostrare graficamente la soluzione di un sistema di tre equazioni di primo grado in tre incognite";

(* Funzioni esercizi *)
exerciseSystemToMatrix::usage = "Mostra l'esercizio in cui viene visualizzato un sistema di equazioni random e l'utente deve trasformarlo nella matrice associata";
exerciseTriangularizeMatrix::usage = "Mostra l'esercizio in cui viene visualizzata una matrice (relativa ad un sistema random) e l'utente deve applicare il Metodo di Gauss per inserire la matrice ridotta associata";
exerciseReducedMatrixToSystem::usage = "Mostra l'esercizio in cui viene visualizzata una matrice ridotta con Gauss e l'utente deve ricostruire il sistema ridotto associato";
exerciseFinalGauss::usage = "Mostra l'esercizio finale in cui l'utente inserisce un sistema di equazioni e applica il Metodo di Gauss passo per passo";
finalExerciseRandomAnswers::usage = "Si occupa di generare le soluzioni random per il sistema dato in input nell'esercizio finale";
finalExerciseFoundAnswer::usage = "Si occupa di concludere l'esercizio finale, ricavando e mostrando la risposta corretta";
diagonalMatrixQuestion::usage = "Domanda: Trovare la diagonale di una matrice";
rowColumnQuestion::usage = "Domanda: Trovare quante righe/colonne possiede la matrice indicata";
indexQuestion::usage = "Domanda: Indentificare un elemento di una matrice, dato l'indice o il valore. ";
identityMatrixQuestion::usage = "Domanda: Viene richiesto dove la matrice identit\[AGrave] presenta degli elementi diversi da zero(Diagonale)";
haveDiagonalQuestion::usage = "Domanda: Viene richiesto quale tra le matrice fornite possiede una diagonale";
questionsExercise::usage = "Gestisce la casualit\[AGrave] delle domande, la correttezza della risposta fornita e la relativa stampa sul Notebook";

(* Funzioni d'appoggio *)
displayEquationSystem::usage = "Permette la visualizzazione di una lista di equazioni sotto forma di sistema";
highlightMatrixElements::usage = "Applica gli stili ad una matrice per evidenziare coefficienti e termini noti";
calculateMatrixDims::usage = "Restituisce la dimensioni (righe, colonne) di una matrice";
transformToMatrix::usage = "Partendo da un sistema di equazioni, restituisce la matrice associata";
getRandomSystem::usage = "Restituisce un sistema di equazioni random tra quelli presenti nel dataset";
oneElementList::usage = "Funzione d'appoggio necessaria per mantenere l'invariante di una lista da un singolo elemento";
stringInputToSystem::usage = "Funzione che si occupa di trasformare un insieme di equazioni in input in un sistema di equazioni valutabili da Mathematica";
checkInputFormat::usage = "Funzione che si occupa di controllare il formato del sistema inserito nell'esercizio finale";

(* variabili riservate ai sistemi *)
t::usage = "";
x::usage = "";
y::usage = "";
z::usage = "";

Begin["`Private`"];

(* Disattivazione warnings *)
Off[Solve::svars];
Off[General::shdw];
Off[Set::setraw];

(* Impostazione della directory del notebook come radice *)
SetDirectory[NotebookDirectory[]];

(* FUNZIONI PRINCIPALI *)

(* Fattorizzazione LU con pivoting parziale;
@PARAM A = matrice quadrata da fattorizzare
@PARAM b = vettore dei termini noti 
Esempio: usata in exerciseTriangularizeMatrix, exerciseReducedMatrixToSystem, exerciseFinalGauss *)
fattorizzazioneLU[A_,b_] := Module[{L,U,P,bPerm,matriceEdited,n,pivot,candidatePivot,subColonna,pivotIndex,lambda,i,k},
	(* se la matrice non \[EGrave] quadrata oppure il determinate \[EGrave] zero, la fattorizzazione LU non \[EGrave] applicabile *)
	If[Not[SquareMatrixQ[A]],Return[{Null,Null,Null,Null}]];
	If[SameQ[Det[A],0],Return[{Null,Null,Null,Null}]];

	(* inizializza n: dimensione della matrice, P: matrice di permutazione, L: matrice triangolare inferiore *)
	n = Length[A];
	matriceEdited = A;
	P = IdentityMatrix[n];
	L = ArrayReshape[ConstantArray[0,n*n],{n, n}]; (* ArrayReshape trasforma una lista di n^2 elementi in una matrice n*n *)

	For[i=1,i<n,i++,
		(* considera gli elementi sotto la diagonale principale (esso compreso) della colonna i-esima *)
		subColonna = matriceEdited[[All,i]][[i;;]]; 
		(* individua il massimo elemento in valore assoluto *)
		candidatePivot = Max[Abs[subColonna]];
		(* ricava la posizione del possibile pivot RELATIVA alla sotto colonna *)
		pivotIndex = FirstPosition[Abs[matriceEdited[[All,i]]][[i;;]],candidatePivot][[1]];
		(* rende la posizione ASSOLUTA rispetto all'intera matrice *)
		pivotIndex = pivotIndex+i-1;
		(* ricava il pivot *)
		pivot = matriceEdited[[pivotIndex,i]];

		If[Not[Equal[i,pivotIndex]],
			(* se il pivot NON \[EGrave] gi\[AGrave] quello sulla diagonale principale, scambia le righe e aggiorna la matrice di permutazione *)
			(* Cycles prende gli indici delle righe da scambiare e, insieme alla Permute, le scambia in modod ciclico *)
			matriceEdited = Permute[matriceEdited,Cycles[{{i,pivotIndex}}]];
			P = Permute[P, Cycles[{Flatten[{i,pivotIndex}]}]];
		];
		For[k=i+1,k<=n,k++,
			If[Not[Equal[matriceEdited[[k,i]],0]],
				(* se il pivot non \[EGrave] 0, calcola lambda, applica la mossa di Gauss a tutti gli elementi della riga k-esima e aggiorna L *)
				lambda = (matriceEdited[[k,i]]/pivot);
				L[[k,i]] = lambda;
				matriceEdited[[k]] = matriceEdited[[k]]-lambda*matriceEdited[[i]];
			];
		];
		(* se \[EGrave] avvenuto uno scambio di righe, anche gli opportuni elementi di L vanno permutati *)
		If[Not[Equal[i,pivotIndex]]&&i>=2,
			L[[i;;n,1;;i-1]]=Permute[L[[i;;n,1;;i-1]],Cycles[{{1,pivotIndex-i+1}}]];
		];
	];
	(* U \[EGrave] la matrice triangolare superiore ottenuta dalla fattorizzazione *)
	U = matriceEdited;
	(* aggiunge la matrice diagonale ad L *)
	L=L+IdentityMatrix[n];
	(* calcola il vettore b aggiornato applicando la formula L^(-1)*P*b *)
	bPerm = Inverse[L].P.b;
	Return[{L,U,P,bPerm}];
];

(* Funzione che permette di mostrare graficamente la soluzione
di un sistema di due equazioni di primo grado in due incognite;
@PARAM eq1 = prima equazione in x e y
@PARAM eq2 = seconda equazione in x e y
Esempio: usata nella seconda slide dei sistemi lineari *)
plotLinearSystem2[eq1_, eq2_] := Module[{solutions, points},
	(* calcola le soluzioni *)
	solutions = Solve[{eq1,eq2}];
	(* crea i punti relativi alle soluzioni *)
	points = {Blue, PointSize[Large], Point[{x,y}/.solutions]};
	(* crea i plot per le rette relative ad ogni equazione *)
	plotEq1 = Plot[y /. Solve[eq1], {x,-10,10},AspectRatio->1,PlotStyle->Red,PlotLegends->{eq1//TraditionalForm}];
	plotEq2 = Plot[y /. Solve[eq2], {x,-10,10},AspectRatio->1,PlotStyle->Orange,PlotLegends->{eq2//TraditionalForm}];
	(* mostra i plot e i punti delle soluzioni *)
	Show[plotEq1, plotEq2, Graphics[{points}], ImageSize->400,Background->White]
];

(* Funzione che permette di mostrare graficamente la soluzione
di un sistema di tre equazioni di primo grado in tre incognite;
@PARAM eq1 = prima equazione in x, y e z
@PARAM eq2 = seconda equazione in x, y e z
@PARAM eq3 = terza equazione in x, y e z
Esempio: usata nella terza slide dei sistemi lineari *)
plotLinearSystem3[eq1_,eq2_,eq3_] := Module[{sols, points},
	(* calcola le soluzioni *)
	sols = Solve[{eq1,eq2,eq3}];
	(* crea i punti relativi alle soluzioni *)
	points = {Red, PointSize[Large], Point[{x,y,z} /. sols]};
	(* crea i plot per i piani relativi ad ogni equazione *)
	plotEq1 = Plot3D[z /. Solve[eq1], {x,-10,10}, {y,-10,10},AspectRatio->1,Mesh->None,PlotStyle->RGBColor[1.0,1.0,0,0.4],PlotLegends->{eq1//TraditionalForm}];
	plotEq2 = Plot3D[z /. Solve[eq2], {x,-10,10}, {y,-10,10},AspectRatio->1,Mesh->None,PlotStyle->RGBColor[0.55,0.87,0.64,0.4],PlotLegends->{eq2//TraditionalForm}];
	plotEq3 = Plot3D[z /. Solve[eq3], {x,-10,10}, {y,-10,10},AspectRatio->1,Mesh->None,PlotStyle->RGBColor[0.99,0.47,1.0,0.4],PlotLegends->{eq3//TraditionalForm}];
	(* mostra i plot e i punti per le soluzioni *)
	Show[plotEq1,plotEq2,plotEq3,Graphics3D[{points}],ImageSize->400]
];

(* FUNZIONI ESERCIZI *)

(* Mostra l'esercizio in cui viene visualizzato un sistema di equazioni (eventualmente random) e l'utente deve trasformarlo nella matrice associata
@PARAM inputSystem = eventuale sistema di equazioni in input
@PARAM result = indica lo step successivo dell'esercizio finale
- viene applicata la HoldRest per mantenere gli ultimi parametri non valutati e simulare un passaggio per riferimento
Esempio: usata nell'esercizio di trasposizione da sistema a matrice e nella verifica finale di Gauss *)
exerciseSystemToMatrix[inputSystem_:Null,result_:0] := DynamicModule[
	{rowCount, colCount, checkButtonSystemToMatrix, restartButton, gridOptions, system, inputMatrixShown, matrix, inputMatrix,
	okColor=RGBColor[0,1,0,0.4], wrongColor=RGBColor[1,0,0,0.4], shown=False, dialogImage, dialogText, output},

	(* se non si ha un sistema in input, ne viene scelto uno random *)
	If[SameQ[inputSystem,Null], system = getRandomSystem[False], system = inputSystem];
	(* calcola la matrice e la sua dimensione *)
	{rowCount,colCount} = calculateMatrixDims[system];
	matrix = transformToMatrix[system,True];
	matrix = Flatten[matrix];
	inputMatrix = ConstantArray[0,rowCount*colCount];
	(* crea una matrice di input field in modo da averne uno associato ad ogni coefficiente/termine noto del sistema *)
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
	(* Button check: controlla che la matrice inserita corrisponda alla soluzione *)
	checkButtonSystemToMatrix = Button[Style["Verifica!",24],
		shown=False;
		(* se la matrice inserita e quella calcolata coincidono, mostra un messaggio di correttezza ed eventualmente passa al passo successivo *)
		(* altrimenti mostra un messaggio d'errore e indica i valori errati *)
		If[MatchQ[matrix,inputMatrix],
			If[Not[SameQ[inputMatrix, Null]], result = 2];
			dialogImage = Import["images/checkmark.png"];
			dialogText = Style["CORRETTO. Bravo!",20,RGBColor[0.14,0.61,0.14]],
			If[Not[SameQ[inputMatrix, Null]], result = 1];
			dialogImage = Import["images/error.png"];
			dialogText = Style["SBAGLIATO. Riprova!",20,Red];
			shown=True;
		];
		MessageDialog[Column[{dialogImage,dialogText}, Spacings->{2,4}, Alignment->Center,
			Frame->All, FrameStyle->RGBColor[0,0,0,0], ItemSize->Fit]],
		ImageSize->{200,50}];
	(* se l'esercizio viene svolto singolarmente, viene creato il pulsante per ricominciare, il quale ricarica la cella del notebook relativa *)
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
	Grid[output, gridOptions]
];
SetAttributes[exerciseSystemToMatrix,HoldRest];

(* Mostra l'esercizio in cui viene visualizzata una matrice (relativa ad un sistema random) e l'utente deve applicare il Metodo di Gauss per inserire la matrice ridotta associata
@PARAM system = sistema in input
@PARAM finalExercise = indica se la funzione viene eseguita come passo dell'esercizio finale o meno
@PARAM result = indica lo step successivo dell'esercizio finale
@PARAM oldResult = indica lo step precedente dell'esericizio finale
- viene applicata la HoldRest per mantenere gli ultimi parametri non valutati e simulare un passaggio per riferimento
Esempio: usata nell'esercizio da matrice a matrice ridotta e nella verifica finale di Gauss *)
exerciseTriangularizeMatrix[system_, finalExercise_:False, result_:0, oldResult_:0] := DynamicModule[
	{systemMatrix, rows, cols, A, b, L, U, P, newB, UFlattened, shown=False, okColor=RGBColor[0,1,0,0.4],
	wrongColor=RGBColor[1,0,0,0.4], inputMatrix, inputMatrixShown, dialogImage, dialogText, checkButtonTriangular, restartButton, gridOptions, output},

	(* ricava la matrice completa associata al sistema e le sue dimensioni *)
	systemMatrix = transformToMatrix[system,True];
	{rows,cols} = calculateMatrixDims[system];
	(* ricava la matrice dei coefficienti attraverso la Drop che toglie l'ultima colonna *)
	A = Drop[systemMatrix, None, -1];
	(* ricava il vettore dei termini noti attraverso la Take che prende l'ultima colonna *)
	b = Take[systemMatrix, All, -1];
	(* applica la fattorizzazione LU su A e b *)
	{L,U,P,newB} = fattorizzazioneLU[A,b];

	gridOptions = {
		Frame->All,
		FrameStyle->RGBColor[0.94,0.94,0.94],
		ItemStyle->Directive[FontFamily->"Roboto Condensed", FontSize->28],
		Spacings->{10,3},
		ItemSize->Fit,
		Alignment->{{Right,Left}}
	};

	(* se si \[EGrave] nell'esercizio finale e la fattorizzazione fallisce, va al passo finale perch\[EGrave] il sistema \[EGrave] impossibile o indeterminato *)
	If[MatchQ[{L,U,P,newB}, ConstantArray[Null,4]], If[finalExercise, result=5],
		(* altrimenti ricava la matrice ridotta completa *)
		UFlattened = Flatten[Join[U,ArrayReshape[newB,{rows,1}],2]];
		inputMatrix = ConstantArray[0, rows*cols];
		(* crea una matrice di input field avente le stesse dimensioni di quella ricavata, in modo da associare ogni input field ad un coefficiente/termine noto *)
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
		(* sostituisce gli input field con degli 0 in tutti gli elementi al di sotto della diagonale principale *)
		For[i=2,i<=Length[inputMatrixShown],i++,
			inputMatrixShown[[i]][[1;;i-1]]=0;
		];
		(* crea il pulsante per la verifica *)
		checkButtonTriangular = Button[Style["Verifica!",24],
			(* shown viene usato per mostrare il colore di sfondo agli input field, in base alla correttezza del contenuto *)
			shown=False;
			(* se la matrice inserita a mano coincide con quella calcolata dalla fattorizzazione, mostra un messaggio di correttezza ed eventualmente passa al passo successivo *)
			(* altrimenti mostra un messaggio di errore e visualizza gli elementi errati *)
			If[MatchQ[UFlattened,inputMatrix],
				If[finalExercise, oldResult = result;result = 3];
				dialogImage = Import["images/checkmark.png"];
				dialogText = Style["CORRETTO. Bravo!",20,RGBColor[0.14,0.61,0.14]],
				If[finalExercise, oldResult = result;result = 2];
				dialogImage = Import["images/error.png"];
				dialogText = Style["SBAGLIATO. Riprova!",20,Red];
				shown=True;
			];
			MessageDialog[Column[{dialogImage,dialogText}, Spacings->{2,4}, Alignment->Center,
				Frame->All, FrameStyle->RGBColor[0,0,0,0], ItemSize->Fit]],
			ImageSize->{200,50}];
			(* crea il pulsante per ricominciare l'esericio se esso viene svolto singolarmente *)
		If[Not[finalExercise],
			restartButton = Button[Style["Ricomincia!",24],
				shown=False;
				NotebookFind[EvaluationNotebook[], "exerciseTriangolarizzazioneCellTag",All,CellTags];
				SelectionEvaluate[EvaluationNotebook[]],
				ImageSize->{200,50}];
		];

		output = {
			{highlightMatrixElements[systemMatrix]//MatrixForm,inputMatrixShown//MatrixForm},
			{checkButtonTriangular,If[Not[finalExercise],restartButton,SpanFromLeft]}
		};
		Grid[output, gridOptions]
	]
];
SetAttributes[exerciseTriangularizeMatrix, HoldRest];

(* Mostra l'esercizio in cui viene visualizzata una matrice ridotta con Gauss e l'utente deve ricostruire il sistema ridotto associato
@PARAM system = eventuale sistema in input
@PARAM result = indica lo step successivo dell'esercizio finale
@PARAM oldResult = indica lo step precedente dell'esericizio finale
- viene applicata la HoldRest per mantenere gli ultimi parametri non valutati e simulare un passaggio per riferimento
Esempio: usata nell'esercizio di scrittura delle equazioni a partire dalla matrice ridotta e nella verifica finale di Gauss *)
exerciseReducedMatrixToSystem[system_:Null,result_:0, oldResult_:0] := DynamicModule[
	{reducedMatrix, inputSystem, inputSystemShown, randomSystem, A, b, U, L, P, newB, matrix, rows, checkButtonReduced, cols, systemToCheck, output,
	okColor=RGBColor[0,1,0,0.4], wrongColor = RGBColor[1,0,0,0.4], shown = False,matrixToCheck, restartButton, gridOptions, dialogImage, dialogText},

	(* se non c'\[EGrave] un sistema in input, ne sceglie uno random; poi calcola la matrice associata e le sue dimensioni *)
	If[SameQ[system, Null], randomSystem = getRandomSystem[], randomSystem = system];
	{rows, cols} = calculateMatrixDims[randomSystem];
	matrix = transformToMatrix[randomSystem,True];

	(* ricava la matrice A ed il vettore dei termini noti *)
	A = Drop[matrix, None, -1];
	b = Take[matrix, All, -1];
	(* applica la fattorizzazione LU su A e b *)
	{L,U,P,newB} = fattorizzazioneLU[A,b];
	(* ottiene la matrice ridotta completa unendo la matrice U e il nuovo vettore b *)
	reducedMatrix = Join[U, ArrayReshape[newB,{rows,1}],2];
	systemToCheck = ConstantArray[Null, rows];
	matrixToCheck = ArrayReshape[ConstantArray[0,rows*cols],{rows, cols}];

	inputSystem = ConstantArray[Null, rows];
	(* crea gli input field, uno per ogni equazione da inserire, che cambiano colore in base alla correttezza del contenuto *)
	inputSystemShown = Partition[
			Table[With[{i=i},
				Dynamic[Framed[
					InputField[Dynamic[inputSystem[[i]]], String,FieldSize->{15,1},Appearance->"Frameless",
						Background->If[shown,
							If[SameQ[reducedMatrix[[i]],matrixToCheck[[i]]],RGBColor[0,1,0,0.4],RGBColor[1,0,0,0.4]],
							White]
						],
						FrameStyle->If[shown,
							If[SameQ[reducedMatrix[[i]],matrixToCheck[[i]]],RGBColor[0,1,0,0.4],RGBColor[1,0,0,0.4]],
							Automatic],
						FrameMargins->None
					]]
				], {i,rows}
			],1];
		(* crea il pulsante di verifica *)
	checkButtonReduced = Button[Style["Verifica!",24],
		If[Not[MemberQ[inputSystem, Null]],
			(* se sono state inserite tutte le equazioni, manipola le stringhe per trasformarle in espressioni valutabili *)
			systemToCheck = StringReplace[#, "=" -> "\[Equal]"] & @ inputSystem;
			systemToCheck = ToExpression[#] & @ systemToCheck;
			(* ricava la matrice associata alle equazioni scritte a mano *)
			matrixToCheck = transformToMatrix[systemToCheck, True];
			(* se la matrice calcolata e quella ricavata dalle equazioni scritte a mano coincidono, mostra il messaggio di correttezza ed eventualmente passa allo step successivo *)
			(* altrimenti mostra un messaggio di errore *)
			If[MatchQ[reducedMatrix,matrixToCheck],
				If[Not[SameQ[system, Null]], oldResult = result;result = 4];
				dialogImage = Import["images/checkmark.png"];
				dialogText = Style["CORRETTO. Bravo!",20,RGBColor[0.14,0.61,0.14]],
				If[Not[SameQ[system, Null]], oldResult=result;result = 3];
				dialogImage = Import["images/error.png"];
				dialogText = Style["SBAGLIATO. Riprova!",20,Red];
			],
			(* mostra un messaggio di errore se non sono state scritte tutte le equazioni *)
			If[Not[SameQ[system, Null]], oldResult=result;result = 3];
			dialogImage = Import["images/error.png"];
			dialogText = Style["SBAGLIATO. Riprova!",20,Red];
		];
		MessageDialog[Column[{dialogImage,dialogText}, Spacings->{2,4}, Alignment->Center,
				Frame->All, FrameStyle->RGBColor[0,0,0,0], ItemSize->Fit]];
		shown = True,
		ImageSize->{200,50}
	];
	If[SameQ[system, Null],
		(* crea il pulsante di restart se si sta svolgendo questo esercizio singolarmente *)
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
		{highlightMatrixElements[reducedMatrix]//MatrixForm,displayEquationSystem[Flatten[inputSystemShown]]},
		{checkButtonReduced,If[SameQ[system, Null],restartButton,SpanFromLeft]}
	};
	Grid[output, gridOptions]
];
SetAttributes[exerciseReducedMatrixToSystem, HoldRest];

(* Funzione che mostra l'esercizio finale in cui l'utente inserisce un sistema di equazioni e applica il Metodo di Gauss passo per passo
Esempio: usata nella slide della verifica finale di Gauss *)
exerciseFinalGauss[] := DynamicModule[
	{inputList = {}, startButton, restartButton, step=0, out, oldStep=0},

	(* crea i pulsanti per iniziare l'esercizio o per ricominciare *)
	startButton = Dynamic[Button[Style["Inizia!",20],
		If[checkInputFormat[inputList[[1]]],
			stringInputToSystem[inputList]; step=1,
			MessageDialog[Column[{Import["images/error.png"],Style["C'\[EGrave] un errore nel sistema inserito.\nPuoi scrivere lettere, numeri, spazi, le operazioni fondamentati (+, -, *, /), l'uguale e la virgola.", 18, TextAlignment->Left]},
				Spacings->{2,4}, Alignment->Center, Frame->All, FrameStyle->RGBColor[0,0,0,0], ItemSize->Fit]];step=0
		], Enabled->If[step!=0,False,True],ImageSize->{180,40}]];
	restartButton = Button[Style["Ricomincia",20], inputList = {}; step=0,ImageSize->{180,40}];

	(* mostra la grid dinamica che cambia contenuto in base al passo corrente dell'esercizio (in base all'opportuna chiamata fatta nello switch);
	inoltre, dopo il primo passo, l'input viene disabilitato per prevenire interferenze *)
	Grid[{
		{Dynamic[InputField[Dynamic[Null, oneElementList[inputList,#]&], String, FieldSize->{40,2}, Enabled->If[step!=0, False,True]]],SpanFromLeft},
		{Dynamic[startButton],restartButton,SpanFromLeft},
		{Dynamic[
			Switch[step,
				0, "",
				1, "Passo 1: Trasforma il sistema nella matrice associata",
				2, "Passo 2: Applica il Metodo di Gauss per ridurre la matrice a gradini",
				3, "Passo 3: Ricava il sistema associato alla matrice ridotta",
				4, "Passo 4: Scegli la soluzione corretta",
				5, "Sei arrivato alla fine!"]
		],SpanFromLeft},
		{Dynamic[
			If[step!=oldStep,
			Switch[step,
				0, "",
				1, exerciseSystemToMatrix[inputList, step],
				2, exerciseTriangularizeMatrix[inputList, True, step,oldStep],
				3, exerciseReducedMatrixToSystem[inputList, step,oldStep],
				4, finalExerciseRandomAnswers[inputList, step,oldStep],
				5, finalExerciseFoundAnswer[inputList]],""]
		,TrackedSymbols:>{step}],SpanFromLeft}
	},
	Frame->All, FrameStyle->RGBColor[0.94, 0.94, 0.94], Spacings->{0,2}, ItemStyle->Directive[FontFamily->"Roboto Condensed", FontSize->28],
	ItemSize->{{Scaled[0.5],Scaled[0.5]}},Alignment->Center]
];

(* Funzione che si occupa di generare le soluzioni random per il sistema dato in input nell'esercizio finale;
@PARAM system = sistema di equazioni dato in input
@PARAM result = indica lo step successivo dell'esercizio finale
@PARAM oldResult = indica lo step precedente dell'esericizio finale
- viene applicata la HoldRest per mantenere gli ultimi parametri non valutati e simulare un passaggio per riferimento
Esempio: usata nella verifica finale di Gauss *)
finalExerciseRandomAnswers[system_, result_:0, oldResult_:0] := DynamicModule[
	{solutions, lhss, variables, answer, answers={}, randomAnswer, radioAnswers, choice, checkButton, gridOptions, i, k, dialogImage, dialogText},

	(* calcola la soluzione reale *)
	solutions = Solve[system];
	(* attraverso Level[#, 1][[1]] e la relativa Map (/@), prendo tutte le parti a sinistra dell'uguale di ogni equazione del sistema *)
	(* per poter estrapolarne le incognite con la builtin Variables *)
	lhss = Level[#, 1][[1]]& /@ system;
	variables = Variables[lhss];
	answer = StandardForm[#] == (# /. solutions[[1]])& /@ variables;
	(* genera una lista di 3 risposte random errate *)
	For[i=1, i<=3, i++,
		randomAnswer = {};
		For[k=1, k<=Length[variables], k++,
			(* le risposte random possono essere semplici numeri oppure frazioni, sia positivi che negativi *)
			AppendTo[randomAnswer, StandardForm[variables[[k]]]==RandomChoice[{RandomInteger[{-20,20}]/RandomInteger[{1,20}], RandomInteger[{-20,20}]}]];
		];
		AppendTo[answers, randomAnswer];
	];
	(* permuta a caso le quattro risposte e trasformale in sistema *)
	answers = Permute[AppendTo[answers,answer], RandomPermutation[4]];
	answers = displayEquationSystem[#]& /@ answers;
	(* crea i radio button associati alle risposte *)
	radioAnswers = Dynamic[RadioButtonBar[Dynamic[choice],answers, Appearance->"Horizontal"]];
	(* crea il pulsante per controllare la risposta scelta *)
	checkButton = Button[Style["Verifica!",24],
		(* se la risposta \[EGrave] giusta, l'utente pu\[OGrave] andare al passo successivo dell'esercizio; altrimenti rimane in quello corrente e viene mostrato un errore *)
		If[SameQ[choice,displayEquationSystem[answer]], oldResult = result; result = 5,
			oldResult = result; result = 4;
			dialogImage = Import["images/error.png"];
			dialogText = Style["SBAGLIATO. Riprova!",20,Red];
			MessageDialog[Column[{dialogImage,dialogText}, Spacings->{2,4}, Alignment->Center,
				Frame->All, FrameStyle->RGBColor[0,0,0,0], ItemSize->Fit]]
		],
		ImageSize->{200,50}];
	gridOptions = {
		Frame->All,
		FrameStyle->RGBColor[0.94,0.94,0.94],
		ItemStyle->Directive[FontFamily->"Roboto Condensed", FontSize->24],
		Spacings->{10,3},
		ItemSize->Fit,
		Alignment->Center
	};
	(* mostra le risposte *)
	Grid[{{radioAnswers,SpanFromLeft},{checkButton,SpanFromLeft}},gridOptions]
];
SetAttributes[finalExerciseRandomAnswers, HoldRest];

(* Funzione che si occupa di concludere l'esercizio finale, ricavando e mostrando la risposta corretta;
@PARAM system = sistema di equazioni utilizzate nell'esercizio
Esempio: usata nella verifica finale di Gauss *)
finalExerciseFoundAnswer[system_] := Module[{solutions, lhss, variables, answer, gridOptions,coeffMatrix},
	(* calcola le soluzioni e ricava la matrice dei coefficienti *)
	solutions = Solve[system];
	coeffMatrix = Drop[transformToMatrix[system,True], None, -1];

	gridOptions = {
		Frame->All,
		FrameStyle->RGBColor[0.94,0.94,0.94],
		ItemStyle->Directive[FontFamily->"Roboto Condensed", FontSize->28],
		Spacings->{10,3},
		ItemSize->Fit,
		Alignment->Center
	};

	If[SquareMatrixQ[coeffMatrix]&&Det[coeffMatrix]!=0,
		(* se il sistema \[EGrave] risolvibile, ricava la soluzione e la mostra sotto forma di sistema; altrimenti mostra il sistema in input indicando se \[EGrave] impossibile o indeterminato *)
		lhss = Level[#, 1][[1]]& /@ system;
		variables = Variables[lhss];
		answer = StandardForm[#] == (# /. solutions[[1]])& /@ variables;
		Grid[{{"Hai risolto correttamente l'intero esercizio!\nIl Metodo di Gauss non ha pi\[UGrave] segreti per te!"},{displayEquationSystem[answer]}},gridOptions],
		Grid[{{displayEquationSystem[system], If[SameQ[Solve[system],{}],"Sistema impossibile","Sistema indeterminato"]}},gridOptions]
	]
];

(* Data una matrice generata random viene richiesto all'utente di selezionarne la diagonale
Esempio: usata nella funzione questionsExercise *)
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
	Return[{text,answers,solution,matrix}]
];

(*Funzione che data una matrice di grandezza casuale richiede all'utente quale delle due indica la colonna o la riga
Esempio: usata nella funzione questionsExercise *)
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

(*Domanda: Indentificare un elemento di una matrice, dato l'indice o il valore
Esempio: usata nella funzione questionsExercise *)
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

(*Domanda che richiede all'utente la caratteristia chiave delle matrici identit\[AGrave]
Esempio: usata nella funzione questionsExercise *)
identityMatrixQuestion[]:=Module[{text,solution,wrongAnswers,answers},
	text = "La matrice identit\[AGrave] presenta degli elementi diversi da zero:";
	solution = "Nella diagonale";
	wrongAnswers = {"Nella prima riga","Nella prima colonna","Indifferente, purch\[EGrave] siano 1"};
	(*Vengono unite la soluzione e le risposte errate*)
	answers = Flatten[{solution,wrongAnswers}];
	Return[{text,answers,solution,Null}]
];

(*Domanda che richiede all'utente quale matrice possiede una diagonale
Esempio: usata nella funzione questionsExercise *)
haveDiagonalQuestion[]:= Module[{text,solution,wrongAnswers,answers},
	text = "Quale delle seguenti matrici possiede una DIAGONALE ?";
	solution = "Matrice quadrata";
	wrongAnswers = {"Matrice rettangolare","Entrambe","Nessuna"};
	(*Vengono unite la soluzione e le risposte errate*)
	answers = Flatten[{solution,wrongAnswers}];
	Return[{text,answers,solution,Null}]
];

(*Funzione che sceglie random le domande a risposta multipla e ne gestisce la relativa stampa sul Notebook
Esempio: usata nella slide dell'esercizio sulle domande di teoria *)
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
	(*If[MatchQ[question,diagonalMatrixQuestion],
		appearance="Horizontal",
		appearance="Vertical"
	];*)
	radioAnswers= Dynamic[RadioButtonBar[Dynamic[choice],answers, Appearance->"Vertical"]];
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

(* FUNZIONI D'APPOGGIO *)

(* Funzione che data in input una lista di equazione la visualizza sotto forma di sistema;
@PARAM eqs = lista di equazioni (espressioni) 
Esempio: usata in tutte le slide in cui va mostrato un sistema lineare *)
displayEquationSystem[eqs_] := Module[{eqsFF},
	(* applica prima HoldForm e poi TraditionalForm a ogni equazione *)
	eqsFF = TraditionalForm /@ HoldForm /@ eqs;
	(* mostra il sistema con la classica visualizzazione matematica *)
	DisplayForm@RowBox[{StyleBox["{", SpanMaxSize->Infinity], Column[eqsFF, Alignment->Left]}]
];

(* Funzione che permette di evidenziare gli elementi di una matrice;
@PARAM matrix = matrice da evidenziare
Esempio: usata in tutte le slide in cui si evidenziano i coefficienti ed i termini noti di una matrice completa *)
highlightMatrixElements[matrix_]:= Module[{lastCol,editedMatrix},
	lastCol = 1;
	(* calcola numero massimo di colonne *)
	For[i=1,i<=Length[matrix],i++,
		If[Length[matrix[[i]]] > lastCol,lastCol = Length[matrix[[i]]]];
	];
	(* colora i coefficienti della matrice di blu (dal primo al penultimo elemento di ogni riga) *)
	editedMatrix = MapAt[Style[#,RGBColor[0.13,0.52,0.96]]&,matrix,{{Range[1,Length[matrix]],Range[1,lastCol-1]}}];
	(* colora i termini noti della matrice di verde (l'ultimo elemento per ogni riga) *)
	editedMatrix = MapAt[Style[#,RGBColor[0.14,0.61,0.14]]&,editedMatrix,{{Range[1,Length[matrix]],lastCol}}];
	Return[editedMatrix];
];

(* Funzione che restituisce la dimensioni (righe, colonne) di una matrice, a partire da un sistema;
@PARAM system = sistema di equazioni da cui ricavare la matrice
Esempio: usata nelle funzioni degli esercizi *)
calculateMatrixDims[system_] := Module[{rows,cols,lhsParts},
	rows = Length[system];
	(* per ogni equazione, considera la parte senza termine noto *)
	lhsParts = Level[#,1][[1]] & /@ system;
	(* calcola il numero di variabili, incrementate di 1, per ottenere il numero di colonne *)
	cols = Length[Variables[lhsParts]]+1;
	Return[{rows,cols}];
];

(* Funzione che prende in input un sistema di equazioni e restituisce la matrice associata;
@PARAM eqs = lista di equazioni (espressioni)
@PARAM withTerms = indica se le equazioni hanno termini noti
Esempio: usata nelle funzioni degli esercizi ed in alcune slide in cui si mostra la matrice associata ad un sistema *)
transformToMatrix[eqs_,withTerms_:False] := Module[{system,matrix,terms,row,incognite,rules,i,j,key},
	If[withTerms,
		(* se withTerms = True, si ricavano i termini noti e la parte dell'equazione senza di essi *)
		system = Level[#,1][[1]] & /@ eqs;
		terms = Level[#,1][[2]] & /@ eqs,
		(* altrimenti le espressioni sono gi\[AGrave] senza termini noti *)
		system = eqs
	];
	(* ricava le incognite *)
	incognite = Sort[Variables[system]];
	(* ricava le associazioni tra incognita e coefficiente (nella forma {0,1,0}->3, ad esempio) *)
	rules = CoefficientRules[system,incognite];
	matrix = {};
	(* crea la matrice *)
	For[i=1,i<=Length[rules],i++,
		row = {};
		(* crea ogni riga *)
		For[j=1,j<=Length[incognite],j++,
			(* crea la chiave per ogni incognita, per risalire al suo coefficiente *)
			key = Normal[SparseArray[{j->1},Length[incognite]]];
			(* se la chiave \[EGrave] nella lista delle regole, allora ricava il coefficiente e lo mette nella riga; altrimenti mette 0 *)
			If[MemberQ[Keys[rules[[i]]],key],
				AppendTo[row,Lookup[rules[[i]],Key[key]]],
				AppendTo[row,0]
			];
		];
		(* se presenti, aggiunge ogni termine noto alla fine di ogni riga *)
		If[withTerms,AppendTo[row,terms[[i]]]];
		(* inserisce la riga creata nella matrice *)
		AppendTo[matrix,row];
	];
	Return[matrix];
];

(* Restituisce un sistema random tra quelli presenti. Funzione utilizzata negli esercizi;
@PARAM solvable = permette di discrimiare tra sistemi risolvibili o anche non risolvibili 
Esempio: usata nelle funzioni degli esercizi *)
getRandomSystem[solvable_:True] := Module[{systems},
	(* definisce una lista di sistemi risolvibili *)
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
	(* se solvable = False, aggiunge anche sistemi impossibili/indeterminati *)
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
	(* ritorna un sistema random tra i precedenti *)
	Return[RandomChoice[systems]];
];

(* Funzione che garantisce la presenza di un unico elemento in una lista;
@PARAM list = lista che deve contenere un solo elemento
@PARAM element = elemento da inserire
- viene usato HoldAll per non valutare i parametri e simulare un passaggio per riferimento
Esempio: usata nella funzione dell'esercizio finale di Gauss *)
oneElementList[list_, element_] := Module[{},
	(* se list ha gi\[AGrave] un elemento viene resettata *)
	If[Not[SameQ[list, {}]], list = {}];
	(* il nuovo elemento viene inserito *)
	list = Insert[list, element, 1];
];
SetAttributes[oneElementList, HoldAll];

(* Funzione che manipola una stringa di equazioni separate da virgola e le trasforma in equazioni valutabili
@PARAM list = lista di equazioni da manipolare
- viene usato HoldAll per non valutare il parametro e simulare un passaggio per riferimento
Esempio: usata nella funzione dell'esercizio finale di Gauss *)
stringInputToSystem[list_] := Module[{},
	(* per ogni equazione, trasforma il singolo = in quello doppio necessario per Mathematica *)
	list = StringReplace[#,"="->"\[Equal]"]& @ list;
	(* splitta la stringa nelle singole equazioni e poi le trasforma in espressioni valutabili *)
	list = ToExpression @ Flatten @ StringSplit[#, ","] & @ list;
];
SetAttributes[stringInputToSystem, HoldAll];

(* Funzione che si occupa di controllare il formato del sistema inserito nell'esercizio finale
@PARAM stringInput = sistema di equazioni scritto dall'utente sotto forma di stringa
Esempio: usata nella funzione dell'esercizio finale di Gauss *)
checkInputFormat[stringInput_] := Module[{inputChars, allowedChars},
	inputChars = Characters[stringInput];
	allowedChars = CharacterRange["a", "z"];
	allowedChars = Join[allowedChars, CharacterRange["A", "Z"]];
	allowedChars = Join[allowedChars, CharacterRange["0", "9"]];
	allowedChars = Join[allowedChars, {"+", "-", "*", "/", ",", " ", "="}];
	Return[AllTrue[inputChars, MemberQ[allowedChars, #] &]];
];

End[];
Protect["GaussLinearSystemsPackage`*"]
EndPackage[];
