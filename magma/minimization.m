/* 
* Developped with Andreas Pieper (https://andreas-pieper.github.io/)
*/

function Hessian(f)
	R := Parent(f);
	return Determinant(Matrix([[Derivative(Derivative(f, 1, i), 1, j) : j in [1..Rank(R)]] : i in [1..Rank(R)]]));
end function;

function QuadraticFormToMatrix(Q)
	R := Parent(Q);
	K := CoefficientRing(R);
	Q_mat := Matrix(K, [[MonomialCoefficient(Q, R.i*R.j) : j in [1..4]] : i in [1..4]]);

	for i in [1..4] do
		Q_mat[i][i] := 2*Q_mat[i][i];
	end for;

	return Q_mat;
end function;

function NewBasis(Q)
	Q_mat := QuadraticFormToMatrix(Q);
	K := BaseRing(Parent(Q));
	M1 := KMatrixSpace(K,4,4);
	D, P, t := OrthogonalizeGram(Q_mat);
	if t lt 3 then
		"The quadric is not of rank 3 or 4";
		return D, t;

	elif t eq 4 then
		L := [-D[4][4]/D[1][1], -D[3][3]/D[2][2]];

		if Category(K) ne FldCom and Category(K) ne FldAC then
			S := AlgebraicClosure(K);
		else
			S := K;
		end if;

		M2 := KMatrixSpace(S,4,4);
		P := ChangeRing(P, S);
		P_fin := 2*Sqrt(S!L[1])*Sqrt(S!L[2])*D[1][1]*D[2][2]*(M2![S!1/(2*D[1][1]),0,0,S!1/(2*D[1][1]*Sqrt(S!L[1])),0,S!-1/(2*D[2][2]),S!-1/(2*D[2][2]*Sqrt(S!L[2])),0,0,S!1/2,S!-1/(2*Sqrt(S!L[2])),0,S!1/2,0,0,S!-1/(2*Sqrt(S!L[1]))])*P;
		return P_fin, 4;

	else
		i := 1;
		while D[i][i] ne 0 do
			i := i+1;
		end while;

		L_swap := [1,2,3,4];
		L_swap[i] := 4;
		L_swap[4] := i;
		P_swap := PermutationMatrix(K, L_swap);
		
		D := P_swap*D*P_swap;
		P := P_swap*P;
		L := [-D[3][3]/D[1][1]];
		
		if Category(K) ne FldCom and Category(K) ne FldAC then
			S := AlgebraicClosure(K);
		else
			S := K;
		end if;
		
		M2 := KMatrixSpace(S,4,4);
		P := ChangeRing(P, S);
		P_fin := 2*Sqrt(S!L[1])*(M2![S!-D[2][2]/(2*D[1][1]),0,S!-D[2][2]/(2*D[1][1]*Sqrt(S!L[1])),0,0,1,0,0,S!1/2,0,S!-1/(2*Sqrt(S!L[1])),0,0,0,0,S!1])*P;
		return P_fin, 3;
	end if;
end function;

function ChangeOfBasis(C, P)
	R := ChangeRing(Parent(C), Parent(P[1][1]));
	C := R!C;
	P := Transpose(Matrix([[R!P[i][j] : j in [1..NumberOfColumns(P)]] : i in [1..NumberOfRows(P)]]));
	Y := ElementToSequence(P*Matrix([[R.i] : i in [1..Rank(R)]]));
	return Evaluate(C, Y);
end function;
	
function CubicNewBasis(Q, C)
	R := Parent(C);
	P := NewBasis(Q);
	R := ChangeRing(R, BaseRing(Parent(P)));
	C1 := R!C;
	return ChangeOfBasis(C1, P);
end function;

function Transvectant3(f, g, h, d)
    P := Parent(f);
    if f eq 0 or g eq 0 or h eq 0 then return P!0; end if;
    
    Sf := [[[Derivative(Derivative(Derivative(f, i, 1), j, 2), k, 3) : k in [0..d]] : j in [0..d]] : i in [0..d]];
	Sg := [[[Derivative(Derivative(Derivative(g, i, 1), j, 2), k, 3) : k in [0..d]] : j in [0..d]] : i in [0..d]];
    Sh := [[[Derivative(Derivative(Derivative(h, i, 1), j, 2), k, 3) : k in [0..d]] : j in [0..d]] : i in [0..d]];

    Tfgh := P!0;
    
    for i := 0 to d do
		for j := 0 to d-i do
            for k := 0 to d-i-j do
                for l := 0 to d-i-j-k do
                    for m := 0 to d-i-j-k-l do
                        n := d-i-j-k-l-m;
                        Tfgh +:= (-1)^(l+m+n)*IntegerRing()!(Factorial(d)/(Factorial(i)*Factorial(j)*Factorial(k)*Factorial(l)*Factorial(m)*Factorial(n)))*(Sf[i+l+1][j+m+1][k+n+1]*Sg[k+m+1][i+n+1][j+l+1]*Sh[j+n+1][k+l+1][i+m+1]);
                    end for;
                end for;
            end for;
		end for;
    end for;

	return Tfgh;
end function;

function DPairing(P, Q)
    R0 := Parent(P);
    S0 := Parent(Q);
    Coef, Mon := CoefficientsAndMonomials(P);
    Exp := [[Degree(m, i) : i in [1..Rank(R0)]] : m in Mon];
    res := 0;
    for j in [1..#Exp] do
        f0 := Q;
        for i in [1..Rank(S0)] do
            f0 := Derivative(f0, Exp[j][i], S0.i);
        end for;
        res +:= Coef[j]*f0;
    end for;
    return res;
end function;

function MinimizeReduceQuadric(Q)
	g := Nrows(Q);
	S := IdentityMatrix(Rationals(), g);
	det := Determinant(Q);
	if det ne 0 then
		fac := Factorization(det);
		primes := [fa: fa in fac|fa[2] gt 1];
		for pri in primes do
			p := pri[1];
			ei := pri[2];
			while ei gt 1 do
				Qred := ChangeRing(Q, GF(p));
				ker := KernelMatrix(Qred);
				nro := Nrows(ker);
				if nro ge ei then
					continue pri;
				end if;
				kerZ := ChangeRing(ker, Integers());
				D, U, V := SmithForm(kerZ);
				SNF := V^-1;
				S := SNF * S;
				Q := SNF * Q * Transpose(SNF);
				j := [i: i in [1..nro]| Valuation(Q[i,i],p) gt 1];
				if j ne [] then
					j := j[1];
				else
					Qsub := Submatrix(Q, 1, 1, nro, nro) div p;
					Qsub := ChangeRing(Qsub, GF(p));
					kersub := KernelMatrix(Qsub);
					kersubZ := ChangeRing(kersub, Integers());
								D, U, V := SmithForm(kersubZ);
								SNF := V^-1;
					Stra := Matrix(BlockDiagMat(<SNF, IdentityMatrix(Integers(), g-nro)>));
					S := Stra * S;
					Q := Stra * Q * Transpose(Stra);
					j := 1;
				end if;
				Dia := DiagonalMatrix([1: i in [1..j-1]] cat [1/p] cat [1: i in [1..g-j]]);
							S := Dia * S;
							Q := Dia * ChangeRing(Q, Rationals()) * Dia;
							Q := ChangeRing(Q, Integers());
							ei -:= 2;
			end while;
		end for;
		R := PolynomialRing(Rationals(), g);
		X := Matrix([[R.i: i in [1..g]]]);
		_, S1 := LLLGram(Q);
		Q := S1 * Q * Transpose(S1);
		return Q, S1*S;
	else
		// A CODER !!
		return Q, IdentityMatrix(Rationals(), 4);
	end if;
end function;

intrinsic MinimizeG4(Q::RngMPolElt, Gamma::RngMPolElt) -> RngMPolElt, RngMPolElt
{Given a non-hyperelliptic genus 4 curve given as the locus of a homogeneous quadratic form Q and a homogeneous cubic form Gamma in 4 variables over Q or Z, returns a model with smaller coefficients}
	P3 := ProjectiveSpace(Rationals(), 3);
	R4<X,Y,Z,T> := CoordinateRing(P3);
	Q := R4!Evaluate(Q, [X,Y,Z,T]);
	Gamma := R4!Evaluate(Gamma, [X,Y,Z,T]);
	Qmat := ZeroMatrix(Rationals(), 4);
	for i in [1..4] do
		for j in [i..4] do
			Qmat[i,j] := MonomialCoefficient(Q, R4.i*R4.j);
		end for;
	end for;
	Qmat := Qmat+Transpose(Qmat);
	Qmat *:= LCM([Denominator(q): q in Eltseq(Qmat)]);
        Qmat /:= GCD([Numerator(e) : e in Eltseq(Qmat)]);
	Qmat := ChangeRing(Qmat, Integers());
	Qmat1, S := MinimizeReduceQuadric(Qmat);
        L := Matrix([[R4.i: i in [1..4]]]);
	Gamma1 := Evaluate(Gamma, Eltseq(L*ChangeRing(S, R4)));
	Gamma1 *:= LCM([Denominator(q): q in Coefficients(Gamma1)]);
        Gamma1 /:= GCD([Numerator(q): q in Coefficients(Gamma1)]);
	Q1 := (L * ChangeRing(Qmat1, R4) * Transpose(L))[1,1];
        Q1 /:= GCD([Numerator(q): q in Coefficients(Q1)]);

	f0 := CubicNewBasis(Q1,Gamma1);
	if Rank(QuadraticFormToMatrix(Q)) eq 3 then
		R<s, t, w> := PolynomialRing(BaseRing(Parent(f0)), [1,1,2]);
		f_weighted := Evaluate(f0, [s^2, s*t, t^2, w]);

		alpha := MonomialCoefficient(f_weighted, w^3);
		f_weighted /:= alpha;        
		f_weighted := Evaluate(f_weighted, [s, t, w-ExactQuotient(Terms(f_weighted, w)[3], 3*w^2)]);

		r1 := BaseRing(Parent(f_weighted)).1;
		f_weighted := Evaluate(f_weighted, [s/r1, t, w]);
		S<[x]> := PolynomialRing(BaseRing(Parent(f_weighted)), 2);

		f := S!Evaluate(f_weighted, [x[1], x[2], 0]);
		v := S!Evaluate(ExactQuotient(Terms(f_weighted, w)[2], w), [x[1], x[2], 0]);

		h24 := Transvectant(f, f, 4);
		inv := Rationals()!(Evaluate(Transvectant(h24, v, 4), [0,0]));
	else 
		R<x, y, u, v> := PolynomialRing(BaseRing(Parent(f0)), 4);
		f_bic := Evaluate(f0, [x*u, y*u, x*v, y*v]);
		
		// Computation of an invariant of the bicubic form
		//Inv := InvariantsGenus4CurvesRank4(f_bic);
		//inv_deg2 := Numerator(Rationals()!DiscriminantFromInvariantsGenus4(Inv));
		inv := Rationals()!(Transvectant(f_bic, f_bic, 3, 3 : invariant := true));
	end if;
		if inv ne 0 then
			// Minimization of the cubic
			l_p := [l[1] : l in Factorization(Numerator(inv)) cat Factorization(Denominator(inv))];
		//inv;
		//l_p;

		for p in l_p do
			Gamma1, S := MinimizeCubicSurface(Gamma1, p);
			Q1 := Evaluate(Q1, Eltseq(L*Transpose(ChangeRing(S, R4))));
		end for;
	end if;

	Q1 := R4!Q1;
	Gamma1 := R4!Gamma1;

	Rt := PolynomialRing(R4, 4); 

	Qt := Evaluate(Q1, [Rt.1, Rt.2, Rt.3, -Rt.1*R4.1-Rt.2*R4.2-Rt.3*R4.3]);
	Qt := Homogenization(Evaluate(Transvectant3(Qt, Qt, Qt, 2), [0,0,0,0]), T);
	Qt /:= GCD([Numerator(c) : c in Coefficients(Qt)]);

	H := DPairing(Qt, Hessian(Gamma1));

	// Reduction of Q and C using the covariant H
	Clust := Cluster(P3, [H, Q1, Gamma1]);
	_, Mat := ReduceCluster(Clust : eps := 10^-100);

	Q1 := R4!ChangeOfBasis(Q1, Transpose(Mat));
	Gamma1 := R4!ChangeOfBasis(Gamma1, Transpose(Mat));
	
	Q1 /:= GCD([Numerator(c) : c in Coefficients(Q1)]);
	Gamma1 /:= GCD([Numerator(c) : c in Coefficients(Gamma1)]);

	mons3 := MonomialsOfDegree(R4, 3);
	LMat := Matrix([[MonomialCoefficient(pos, m): m in mons3] : pos in [Gamma1] cat [Q1*R4.i: i in [1..4]]]);
	A, B := LLL(LMat);
	j := Max([j : j in [1..5] | B[j][1] ne 0]);
	vecopt := [A[j][i] : i in [1..20]];	

	Gamma1 := &+[vecopt[i]*mons3[i] : i in [1..20]];

	return Q1, Gamma1;
end intrinsic;
