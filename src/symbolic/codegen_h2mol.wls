
wavefn[r_, \[Theta]_, \[Phi]_] = (1+c*r)*Exp[-\[Alpha]*r];

(* Wave function VB *)
(* wavefnVB[x_, y_, z_] = wavefn[x, y, z, R]*wavefn[x, y, z, -R]; *)

(* Wave function MO*)
(* wavefnVB[x_, y_, z_] = wavefn[x, y, z, R]+wavefn[x, y, z, -R]; *)

(* gradwavefnVB = Grad[wavefn[x, y, z], {x, y, z}]; *)
(* laplacianwavefnVB = Laplacian[wavefn[x, y, z], {x, y, z}]; *)

gradwavefn = Grad[wavefn[r, \[Theta], \[Phi]], {r, \[Theta], \[Phi]}, "Spherical"];
laplacianwavefn = Laplacian[wavefn[r, \[Theta], \[Phi]], {r, \[Theta], \[Phi]}, "Spherical"];

Print[ToString[TeXForm[FullSimplify[gradwavefn]]]];
Print[ToString[TeXForm[FullSimplify[laplacianwavefn]]]];