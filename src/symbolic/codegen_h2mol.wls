
wavefn[x_, y_, z_, R_] = (1+c*Sqrt[(x-R)^2+y^2+z^2])*Exp[-\[Alpha]*Sqrt[(x-R)^2+y^2+z^2]];

(* Wave function VB *)
wavefnVB[x_, y_, z_] = wavefn[x, y, z, R]*wavefn[x, y, z, -R];

(* Wave function MO*)
wavefnVB[x_, y_, z_] = wavefn[x, y, z, R]+wavefn[x, y, z, -R];

gradwavefnVB = Grad[wavefnVB[x, y, z], {x, y, z}];
laplacianwavefnVB = Laplacian[wavefnVB[x, y, z], {x, y, z}];
Print[ToString[TeXForm[FullSimplify[gradwavefnVB]]]];
Print[ToString[TeXForm[FullSimplify[laplacianwavefnVB]]]];