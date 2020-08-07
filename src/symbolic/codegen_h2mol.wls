
wavefn[r_, \[Theta]_, \[Phi]_, R_] = (1+c*r)*Exp[-\[Alpha]*(r-R)];

(* Wave function VB *)
wavefnVB[r_, \[Theta]_, \[Phi]_] = wavefn[r, \[Theta]_, \[Phi]_, R1]*wavefn[r, \[Theta]_, \[Phi]_, R2];

(* Wave function MO*)
wavefnMO[r_, \[Theta]_, \[Phi]_] = wavefn[r, \[Theta]_, \[Phi]_, R1] + wavefn[r, \[Theta]_, \[Phi]_, R2];

gradwavefnVB = Grad[wavefnVB[r, \[Theta], \[Phi]], {r, \[Theta], \[Phi]}, "Spherical"];
laplacianwavefnVB = Laplacian[wavefnVB[r, \[Theta], \[Phi]], {r, \[Theta], \[Phi]}, "Spherical"];
Print[ToString[TeXForm[FullSimplify[gradwavefnVB]]]];
Print[ToString[TeXForm[FullSimplify[laplacianwavefnVB]]]];