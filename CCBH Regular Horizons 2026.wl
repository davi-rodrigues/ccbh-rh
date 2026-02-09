(* ::Package:: *)

(* ::Title:: *)
(*Cosmologically Coupled Black Holes with Regular Horizons:*)
(* field equations derivation*)


(* ::Text:: *)
(*Here we verify and detail the field equations derivation from 2601.03296 [gr-qc] , by Cadoni, Lima, Pitzalis, Rodrigues and Sanna.*)
(**)
(*This code is based on the FTeV package , which can be downloaded from https://github.com/davi-rodrigues/FTeV*)


(* ::Text:: *)
(**)


<< FTeV`

noArg[x_] := ReplaceAll[x, {y_[\[Eta]] -> y, y_[\[Eta], r] -> y}];

$Coordinates = {\[Eta], r, \[Theta], \[Phi]};

$Metric = a[\[Eta]]^2 DiagonalMatrix[{- Exp @ \[Alpha][\[Eta],r], Exp @ \[Beta][\[Eta],r], r^2, r^2 Sin[\[Theta]]^2}];

% // MatrixForm

$Assumptions = {a[\[Eta]] > 0, r > 0, Exp[\[Beta][\[Eta],r]] > 0, Exp[\[Alpha][\[Eta],r]] > 0};


(* ::Section:: *)
(*Christofell symbol, Ricci tensor and Ricci scalar: for cross-check. *)
(*This section can be skipped.*)


tev["Chr"] //noArg // TensorPrint


tev["Ricci"] //noArg // TensorPrint 


tev["RicciS"] // Simplify // noArg


(* ::Section:: *)
(*Einstein equations*)


(* ::Text:: *)
(*In the following, we will derive the field equations.*)


Gdd[] = Simplify[tev["G"]]; (*Einstien tensor with covariant indices*)
G[] = indices[Gdd[], "dd", "ud"] // Simplify; (*Einstein tensor with mixed indices*)

\[ScriptCapitalT][] = DiagonalMatrix[{-\[Rho][\[Eta],r], pp[\[Eta],r], pr[\[Eta],r], pr[\[Eta],r]}]; (*Energy momentum tensor with mixed indices*)
\[ScriptCapitalT][a_,b_] := \[ScriptCapitalT][][[a+1, b+1]];
G[a_,b_] := G[][[a+1, b+1]];

fieldEqs[] = G[] - k \[ScriptCapitalT][];
fieldEqs[a_,b_] := fieldEqs[][[a+1, b+1]];

k = 8\[Pi];

fieldEqs[] // TensorPrint


(*3a*)
fieldEqs[0,1]==0  // Simplify // MultiplySides[#, 1/(a[\[Eta]] r)] & // Expand // noArg

ruleBetaDot = Solve[fieldEqs[0,1]==0 , Derivative[1,0][\[Beta]][\[Eta],r]][[1,1]]; (*rule for \[Beta] replacement. It will be relevant latter*)
replaceBetaDot[x_] := ReplaceAll[x, ruleBetaDot];


(*3b*)
fieldEqs[0,0]== 0  // Simplify // MultiplySides[#,  a[\[Eta]]^2 Exp[\[Alpha][\[Eta],r]]] & // Expand // noArg


(*3c*)
fieldEqs[1,1] == 0 // Simplify // MultiplySides[#,  a[\[Eta]]^2 Exp[\[Beta][\[Eta],r]]] & // Expand //noArg


(* ::Text:: *)
(*For the remainig equations, we use the energy momentum tensor conservation. *)


d\[ScriptCapitalT][] = indices[
	DCov[indices[\[ScriptCapitalT][], "ud", "dd"]][], 
	"ddd", "udd"
] // Simplify; (*Defines the tensor  \[Del]^aSubscript[Subscript[T, b], c]. Note: DCoV is only defined for fully covariant tensors, hence "indices" is used before DCoV here.*)
d\[ScriptCapitalT][a_,b_, c_] := d\[ScriptCapitalT][][[a+1, b+1, c+1]]; 
div\[ScriptCapitalT][] = ToTensor[Contract[d\[ScriptCapitalT][aa, aa, b], aa], b]; (*Contracts the covariant derivative, yielding \[Del]^aSubscript[Subscript[T, a], b]. To avoid overlap with a[\[Eta]], do not contract an index named "a", "aa" works fine. The function Contract is not fully localized.*)
div\[ScriptCapitalT][a_] := div\[ScriptCapitalT][][[a+1]]; 
div\[ScriptCapitalT][] // TensorPrint


(*eq. 3e*)
div\[ScriptCapitalT][1] == 0 // Simplify // MultiplySides[#, 1/(2r)] & // Expand // noArg

rulePr = Solve[div\[ScriptCapitalT][1] == 0, pr[\[Eta],r]][[1,1]]; (*It will be useful latter*)
replacePr[x_] := ReplaceAll[x, rulePr];


(*eq. 3d*)
ruleAlphaPrime = Solve[div\[ScriptCapitalT][1] == 0 // Simplify, Derivative[0,1][\[Alpha]][\[Eta],r]][[1,1]];
replaceAlphaPrime[x_] := ReplaceAll[x, ruleAlphaPrime];

div\[ScriptCapitalT][0] == 0  // replaceBetaDot // replaceAlphaPrime // Simplify // MultiplySides[#, 1/a[\[Eta]]] & // Expand // noArg

