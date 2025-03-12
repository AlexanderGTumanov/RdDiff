(* ::Package:: *)

version=1.1;
Print["RdDiff ",ToString@NumberForm[version,{20,1}]]
Print["by Alexander G Tumanov"]

Clear[x, xx, \[Delta], d]
Attributes[xx] = Orderless;
Attributes[\[Delta]] = Orderless;

Clear[DD]
DD[f_, II_] :=
	Block[{CC=Union@Join[Cases[{f},xx[II[[1]],__],\[Infinity]],Cases[{f},xx[__,II[[1]]],\[Infinity]]]},
		Sum[2D[f,CC[[ii]]]((x@@II)-x[Complement[List@@CC[[ii]],{II[[1]]}][[1]],II[[2]]]),{ii,Length[CC]}]
	] +
	Block[{CC=Union@Cases[{f},x[II[[1]],__],\[Infinity]]},
		Sum[\[Delta][II[[2]],Last[List@@CC[[ii]]]] Coefficient[f,CC[[ii]]],{ii,Length[CC]}]
	] + 2 D[f,xx[II[[1]]]](x@@II)
DD::usage = "DD[expr,{a,i}] computes the partial derivative of expr with respect to \!\(\*SuperscriptBox[SubscriptBox[\(x\), \(a\)], \(i\)]\).";
	
Clear[ParallelDD]
ParallelDD[f_, II_] :=
	Block[{CC=Union@Join[Cases[{f},xx[II[[1]],__],\[Infinity]],Cases[{f},xx[__,II[[1]]],\[Infinity]]]},
		ParallelSum[2D[f,CC[[ii]]]((x@@II)-x[Complement[List@@CC[[ii]],{II[[1]]}][[1]],II[[2]]]),{ii,Length[CC]}]
	] +
	Block[{CC=Union@Cases[{f},x[II[[1]],__],\[Infinity]]},
		ParallelSum[\[Delta][II[[2]],Last[List@@CC[[ii]]]] Coefficient[f,CC[[ii]]],{ii,Length[CC]}]
	] + 2 D[f,xx[II[[1]]]](x@@II)
ParallelDD::usage = "ParallelDD[expr,{a,i}] computes the partial derivative of expr with respect to \!\(\*SuperscriptBox[SubscriptBox[\(x\), \(a\)], \(i\)]\). This function uses parallel evaluation.";

Clear[Contract]
Contract[f_,II_] :=
	 Block[{Lx1=Union@Cases[{f},x[_,II[[1]]],\[Infinity]],Lx2=Union@Cases[{f},x[_,II[[2]]],\[Infinity]],L\[Delta]1=Union@Cases[{f}/.\[Delta]@@II->0,\[Delta][II[[1]],_],\[Infinity]],L\[Delta]2=Union@Cases[{f}/.\[Delta]@@II->0,\[Delta][II[[2]],_],\[Infinity]]},
		Sum[If[TrueQ[Lx1[[ii,1]]==Lx2[[jj,1]]],xx[Lx1[[ii,1]]],1/2(xx[Lx1[[ii,1]]]+xx[Lx2[[jj,1]]]-xx[Lx1[[ii,1]],Lx2[[jj,1]]])](f//Coefficient[#,Lx1[[ii]]]&//Coefficient[#,Lx2[[jj]]]&),{ii,Length[Lx1]},{jj,Length[Lx2]}] +
		Sum[x[Lx1[[ii,1]],Complement[List@@(L\[Delta]2[[jj]]),{II[[2]]}][[1]]](f//Coefficient[#,Lx1[[ii]]]&//Coefficient[#,L\[Delta]2[[jj]]]&),{ii,Length[Lx1]},{jj,Length[L\[Delta]2]}] +
		Sum[x[Lx2[[jj,1]],Complement[List@@(L\[Delta]1[[ii]]),{II[[1]]}][[1]]](f//Coefficient[#,L\[Delta]1[[ii]]]&//Coefficient[#,Lx2[[jj]]]&),{ii,Length[L\[Delta]1]},{jj,Length[Lx2]}] +
		Sum[\[Delta][Complement[List@@(L\[Delta]1[[ii]]),{II[[1]]}][[1]],Complement[List@@(L\[Delta]2[[ii]]),{II[[2]]}][[1]]](f//Coefficient[#,L\[Delta]1[[ii]]]&//Coefficient[#,L\[Delta]2[[jj]]]&),{ii,Length[L\[Delta]1]},{jj,Length[L\[Delta]2]}]
	 ] + d Coefficient[f,\[Delta]@@II]
Contract::usage = "Contract[expr,{i,j}] computes the contraction of expr on i and j indexes.";

Clear[ParallelContract]
ParallelContract[f_,II_] :=
	 Block[{Lx1=Union@Cases[{f},x[_,II[[1]]],\[Infinity]],Lx2=Union@Cases[{f},x[_,II[[2]]],\[Infinity]],L\[Delta]1=Union@Cases[{f}/.\[Delta]@@II->0,\[Delta][II[[1]],_],\[Infinity]],L\[Delta]2=Union@Cases[{f}/.\[Delta]@@II->0,\[Delta][II[[2]],_],\[Infinity]]},
		ParallelSum[If[TrueQ[Lx1[[ii,1]]==Lx2[[jj,1]]],xx[Lx1[[ii,1]]],1/2(xx[Lx1[[ii,1]]]+xx[Lx2[[jj,1]]]-xx[Lx1[[ii,1]],Lx2[[jj,1]]])](f//Coefficient[#,Lx1[[ii]]]&//Coefficient[#,Lx2[[jj]]]&),{ii,Length[Lx1]},{jj,Length[Lx2]}] +
		ParallelSum[x[Lx1[[ii,1]],Complement[List@@(L\[Delta]2[[jj]]),{II[[2]]}][[1]]](f//Coefficient[#,Lx1[[ii]]]&//Coefficient[#,L\[Delta]2[[jj]]]&),{ii,Length[Lx1]},{jj,Length[L\[Delta]2]}] +
		ParallelSum[x[Lx2[[jj,1]],Complement[List@@(L\[Delta]1[[ii]]),{II[[1]]}][[1]]](f//Coefficient[#,L\[Delta]1[[ii]]]&//Coefficient[#,Lx2[[jj]]]&),{ii,Length[L\[Delta]1]},{jj,Length[Lx2]}] +
		ParallelSum[\[Delta][Complement[List@@(L\[Delta]1[[ii]]),{II[[1]]}][[1]],Complement[List@@(L\[Delta]2[[ii]]),{II[[2]]}][[1]]](f//Coefficient[#,L\[Delta]1[[ii]]]&//Coefficient[#,L\[Delta]2[[jj]]]&),{ii,Length[L\[Delta]1]},{jj,Length[L\[Delta]2]}]
	 ] + d Coefficient[f,\[Delta]@@II]
ParallelContract::usage = "ParallelContract[expr,{i,j}] computes the contraction of expr on i and j indexes. This function uses parallel evaluation.";
