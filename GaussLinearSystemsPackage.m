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

x::usage = "";
y::usage = "";

Begin["`Private`"];

Off[Solve::svars];
Off[General::shdw];
SetDirectory[NotebookDirectory[]];

plotLinearSystem2[eq1_, eq2_] := Module[{solutions, points},
	solutions = Solve[{eq1,eq2},{x,y},Reals];
	points = {Blue, PointSize[Medium], Point[{x,y} /. solutions]};
	plotEq1 = Plot[y /. Solve[eq1], {x,-10,10},AspectRatio->1];
	plotEq2 = Plot[y /. Solve[eq2], {x,-10,10},AspectRatio->1];
	Show[plotEq1, plotEq2, Graphics[{points}]]
];

End[];
Protect["GaussLinearSystemsPackage`*"]
EndPackage[];
