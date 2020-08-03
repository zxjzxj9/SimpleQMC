(* Generate the derivative give the trail wavefunction *)

(* Needs["CCodeGenerator`"] *)
wavefn[r_] = (1+c*r)*Exp[-\[Alpha]*r];
d2wd2rxr = D[r*wavefn[r], {r, 2}]
(* dwdr = D[wavefn[r], {r, 2}] *)
Print["Trial function:"];
Print[ToString[TeXForm[FullSimplify[wavefn[r]]]]];
rho = wavefn[r1] ^2 / wavefn[r2]^2;
Print["Ratios of the rhos:"];
Print[ToString[TeXForm[Simplify[rho]]]];
energy = -d2wd2rxr/(2*r*wavefn[r]) - 1/r;
Print["Energy expressions:"];
Print[ToString[TeXForm[FullSimplify[energy]]]];
(* src = Compile[{{x}}, dwdx];
file = CCodeGenerator[{dwdx}, "wavefn"]; *)
(* c = Compile[ {{x}}, x^2];
file = CCodeGenerate[c, "fun"]*)
