(* ::Package:: *)

version=1.0;
Print["Coordinate differentiation ",ToString@NumberForm[version,{20,1}]]
Print["by Alexander G Tumanov"]

Clear[x, xx, \[Delta], d]
Attributes[xx] = Orderless;
Attributes[\[Delta]] = Orderless;

Clear[DD]
DD[f_,II_] :=
	Block[{CC=Union@Join[Cases[{f},xx[II[[1]],__],\[Infinity]],Cases[{f},xx[__,II[[1]]],\[Infinity]]]},
		Sum[2D[f,CC[[ii]]]((x@@II)-x[First[Complement[List@@CC[[ii]],{II[[1]]}]],II[[2]]]),{ii,Length[CC]}]
	] +
	Block[{CC=Union@Cases[{f},x[II[[1]],__],\[Infinity]]},
		If[CC=={},0,Sum[\[Delta][II[[2]],Last[List@@CC[[ii]]]] Coefficient[f,CC[[ii]]],{ii,Length[CC]}]]
	] + 2 D[f,xx[II[[1]]]](x@@II)
	
Clear[ParallelDD]
ParallelDD[f_,II_] :=
	Block[{CC=Union@Join[Cases[{f},xx[II[[1]],__],\[Infinity]],Cases[{f},xx[__,II[[1]]],\[Infinity]]]},
		ParallelSum[2D[f,CC[[ii]]]((x@@II)-x[First[Complement[List@@CC[[ii]],{II[[1]]}]],II[[2]]]),{ii,Length[CC]}]
	] +
	Block[{CC=Union@Cases[{f},x[II[[1]],__],\[Infinity]]},
		If[CC=={},0,ParallelSum[\[Delta][II[[2]],Last[List@@CC[[ii]]]] Coefficient[f,CC[[ii]]],{ii,Length[CC]}]]
	] + 2 D[f,xx[II[[1]]]](x@@II)

Clear[Contract]
Contract[f_,\[Mu]_] :=
	Block[{CC=If[TrueQ[Head[ExpandAll[f/.\[Delta][\[Mu],_]->0]]==Plus],List@@ExpandAll[f/.\[Delta][\[Mu],_]->0],{f/.\[Delta][\[Mu],_]->0}]/.{0}->{}},
		Sum[Block[{QQ=Sort@Cases[CC[[ii]],x[_,\[Mu]],\[Infinity]]},
			(CC[[ii]]/.x[_,\[Mu]]->1)If[Length[QQ]==1,xx@@(QQ/.x[a_,\[Mu]]:>a),1/2 (xx[First[QQ]/.x[a_,\[Mu]]:>a]+xx[Last[QQ]/.x[a_,\[Mu]]:>a]-xx@@(QQ/.x[a_,\[Mu]]:>a))]
		],{ii,Length[CC]}]
	] +
	Block[{CC=If[TrueQ[Head[ExpandAll[(f/.\[Delta][\[Mu],\[Mu]]->0)-(f/.\[Delta][\[Mu],_]->0)]]==Plus],List@@ExpandAll[(f/.\[Delta][\[Mu],\[Mu]]->0)-(f/.\[Delta][\[Mu],_]->0)],{(f/.\[Delta][\[Mu],\[Mu]]->0)-(f/.\[Delta][\[Mu],_]->0)}]/.{0}->{}},
		Sum[Block[{QQ=Select[Flatten[List@@#&/@Cases[CC[[ii]],\[Delta][\[Mu],_],\[Infinity]]],!TrueQ[#==\[Mu]]&]},
			If[Length[QQ]==2,(\[Delta]@@QQ)(CC[[ii]]/.\[Delta][\[Mu],_]->1),If[FreeQ[CC[[ii]],x[_,\[Mu]]],(\[Delta]@@Join[QQ,QQ])(CC[[ii]]/.\[Delta][\[Mu],_]->1),CC[[ii]]/.{\[Delta][\[Mu],_]->1,x[a_,\[Mu]]:>x[a,First[QQ]]}]]
		],{ii,Length[CC]}]
	] + d Coefficient[f,\[Delta][\[Mu],\[Mu]]]/;!ListQ[\[Mu]]
ParallelContract[f_,\[Mu]\[Mu]_] := First@Nest[{ParallelContract[#[[1]],\[Mu]\[Mu][[#[[2]]]]],#[[2]]+1}&,{f,1},Length[\[Mu]\[Mu]]]/;ListQ[\[Mu]\[Mu]]

Clear[ParallelContract]
ParallelContract[f_,\[Mu]_] :=
	Block[{CC=If[TrueQ[Head[ExpandAll[f/.\[Delta][\[Mu],_]->0]]==Plus],List@@ExpandAll[f/.\[Delta][\[Mu],_]->0],{f/.\[Delta][\[Mu],_]->0}]/.{0}->{}},
		ParallelSum[Block[{QQ=Sort@Cases[CC[[ii]],x[_,\[Mu]],\[Infinity]]},
			(CC[[ii]]/.x[_,\[Mu]]->1)If[Length[QQ]==1,xx@@(QQ/.x[a_,\[Mu]]:>a),1/2 (xx[First[QQ]/.x[a_,\[Mu]]:>a]+xx[Last[QQ]/.x[a_,\[Mu]]:>a]-xx@@(QQ/.x[a_,\[Mu]]:>a))]
		],{ii,Length[CC]}]
	] +
	Block[{CC=If[TrueQ[Head[ExpandAll[(f/.\[Delta][\[Mu],\[Mu]]->0)-(f/.\[Delta][\[Mu],_]->0)]]==Plus],List@@ExpandAll[(f/.\[Delta][\[Mu],\[Mu]]->0)-(f/.\[Delta][\[Mu],_]->0)],{(f/.\[Delta][\[Mu],\[Mu]]->0)-(f/.\[Delta][\[Mu],_]->0)}]/.{0}->{}},
		ParallelSum[Block[{QQ=Select[Flatten[List@@#&/@Cases[CC[[ii]],\[Delta][\[Mu],_],\[Infinity]]],!TrueQ[#==\[Mu]]&]},
			If[Length[QQ]==2,(\[Delta]@@QQ)(CC[[ii]]/.\[Delta][\[Mu],_]->1),If[FreeQ[CC[[ii]],x[_,\[Mu]]],(\[Delta]@@Join[QQ,QQ])(CC[[ii]]/.\[Delta][\[Mu],_]->1),CC[[ii]]/.{\[Delta][\[Mu],_]->1,x[a_,\[Mu]]:>x[a,First[QQ]]}]]
		],{ii,Length[CC]}]
	] + d Coefficient[f,\[Delta][\[Mu],\[Mu]]]/;!ListQ[\[Mu]]
ParallelContract[f_,\[Mu]\[Mu]_] := First@Nest[{ParallelContract[#[[1]],\[Mu]\[Mu][[#[[2]]]]],#[[2]]+1}&,{f,1},Length[\[Mu]\[Mu]]]/;ListQ[\[Mu]\[Mu]]
