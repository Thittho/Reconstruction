import "conic_param.m" : ConicParametrization;

function ReconstructionGenus4Rank4(L : minimize := true)
	R<X,Y,Z,T> := PolynomialRing(Parent(L[1]), 4);
	
	c3c3 := L[11];
	c3c51 := L[12];
	c3c52 := L[13];
	c3c53 := -772/9*L[6]-772/27*L[1]*L[4]+965/96*L[1]^4+193/96*L[3]^2+193/72*L[2]*L[3]-193/16*L[1]^2*L[3]-193/72*L[1]^2*L[2]+416/9*L[1]*L[11]+1955/3*L[12]-512/27*L[13];
	c51c51 := L[14];
	c51c52 := L[15];
	c51c53 := L[16];
	c52c52 := L[17];
	c52c53 := L[18];
	c53c53 := L[19];

	M_3_51_52_53 := Matrix(Parent(L[1]), [[c3c3, c3c51, c3c52, c3c53], [c3c51, c51c51, c51c52, c51c53], [c3c52, c51c52, c52c52, c52c53], [c3c53, c51c53, c52c53, c53c53]]);

	if Determinant(M_3_51_52_53) ne 0 then
		fc3c3c3 := L[8];
		fc3c3c51 := L[26];
		fc3c3c52 := L[27];
		fc3c3c53 := L[28];
		fc3c51c51 := L[34];
		fc3c51c52 := L[35];
		fc3c51c53 := L[36];
		fc3c52c52 := L[37];
		fc3c52c53 := L[38];
		fc3c53c53 := L[39];
		fc51c51c51 := L[45];
		fc51c51c52 := L[46];
		fc51c51c53 := L[47];
		fc51c52c52 := L[48];
		fc51c52c53 := L[49];
		fc51c53c53 := L[50];
		fc52c52c52 := L[51];
		fc52c52c53 := L[52];
		fc52c53c53 := L[53];
		fc53c53c53 := L[54];
		
		Q := c3c3*X^2+2*c3c51*X*Y+2*c3c52*X*Z+2*c3c53*X*T+c51c51*Y^2+2*c51c52*Y*Z+2*c51c53*Y*T+c52c52*Z^2+2*c52c53*Z*T+c53c53*T^2;
		E := fc3c3c3*X^3+3*fc3c3c51*X^2*Y+3*fc3c3c52*X^2*Z+3*fc3c3c53*X^2*T+3*fc3c51c51*X*Y^2+6*fc3c51c52*X*Y*Z+6*fc3c51c53*X*Y*T+3*fc3c52c52*X*Z^2+6*fc3c52c53*X*Z*T+3*fc3c53c53*X*T^2+fc51c51c51*Y^3+3*fc51c51c52*Y^2*Z+3*fc51c51c53*Y^2*T+3*fc51c52c52*Y*Z^2+6*fc51c52c53*Y*Z*T+3*fc51c53c53*Y*T^2+fc52c52c52*Z^3+3*fc52c52c53*Z^2*T+3*fc52c53c53*Z*T^2+fc53c53c53*T^3; 

		return Q, E;
	end if;
	
    "Not implemented yet";
    return 0, 0;
end function;

function ReconstructionGenus4Rank3(I : minimize := true)
    a11 := 1/5*I[1]*I[6] + I[9] + 12/7*I[10];
    a12 := I[13];
    a13 := -I[12];
    a22 := I[3];
    a23 := -I[21];
    a33 := 1/2*(-1/2*I[1]*I[10] + 1/2*I[2]*I[6] - I[8]^2 + 3/2*I[22]);

    v11 := 2/35*I[1]*I[7] + 1/7*I[6]*I[8] + 16/11*I[14] + I[15];
    v12 := 2/7*I[22] + I[23];
    v13 := I[16] - 3/2*I[19];
    v22 := I[24];
    v23 := I[25];
    v33 := 1/2*(1/12*I[1]*I[6]*I[8] - 1/2*I[2]*I[7] - 1/4*I[6]*I[13] + I[8]*I[10]);

    f111 := 2/35*I[1]^2*I[7] - 23/1155*I[1]*I[6]*I[8] + 6/11*I[1]*I[14] +
        2/9*I[1]*I[15] - 16/77*I[2]*I[7] + 6/11*I[6]*I[13] - 1/2*I[8]*I[9] -
        52/77*I[8]*I[10] + 5/3*I[27] + 14/11*I[29] - 2*I[30];
    f112 := 1/7*I[1]^2*I[10] - 1/30*I[1]*I[2]*I[6] + 1/6*I[1]*I[8]^2 - 1/7*I[1]*I[22]
        + 2/9*I[1]*I[23] - I[2]*I[9] - 4/7*I[2]*I[10] + 1/5*I[3]*I[6] + 8/3*I[35] +
        6/7*I[36];
    f113 := 1/2*I[1]*I[16] - 2/3*I[1]*I[19] - 1/12*I[2]*I[11] - 3/14*I[6]*I[21] +
        1/2*I[8]*I[12] + 5/4*I[33] - 3/2*I[41];
    f122 := -1/54*I[1]^3*I[8] + 1/18*I[1]^2*I[13] + 1/9*I[1]*I[2]*I[8] +
        1/6*I[1]*I[24] - 1/3*I[2]*I[13] - 1/3*I[3]*I[8] + 2*I[43];
    f123 := -1/18*I[1]^2*I[12] + 10/21*I[1]*I[25] + 1/6*I[1]*I[26] + 1/3*I[2]*I[12] +
        1/3*I[8]*I[21] - I[42];
    f133 := -19/6160*I[1]^2*I[6]*I[8] - 27/220*I[1]*I[2]*I[7] +
        57/6160*I[1]*I[6]*I[13] - 1/8*I[1]*I[8]*I[9] + 123/308*I[1]*I[8]*I[10] -
        3/308*I[2]*I[6]*I[8] + 3/44*I[2]*I[14] + 5/24*I[2]*I[15] + 24/77*I[3]*I[7] -
        3/112*I[6]*I[24] + 3/8*I[8]^3 - 295/308*I[8]*I[22] - 1/4*I[8]*I[23] +
        3/8*I[9]*I[13] - 38/77*I[10]*I[13] + 6/11*I[44] - 5/4*I[48];
    f222 := I[4];
    f223 := -1/12*I[1]^2*I[21] + 1/4*I[1]*I[37] + 1/2*I[2]*I[21] - I[50];
    f233 := 1/72*I[1]^2*I[2]*I[6] - 1/18*I[1]^2*I[8]^2 - 1/24*I[1]*I[3]*I[6] + 1/4*I[1]*I[8]*I[13] - 1/12*I[2]^2*I[6] +
        1/3*I[2]*I[8]^2 - 3/4*I[2]*I[22] + 1/2*I[3]*I[10] - 1/4*I[13]^2;
    f333 := 1/2*(1/6*I[1]*I[2]*I[11] + 1/24*I[1]*I[6]*I[21] - 1/2*I[1]*I[8]*I[12] + 3/2*I[2]*I[16] + 5/4*I[2]*I[19] -
        1/2*I[3]*I[11] - 1/8*I[6]*I[37] - 9/14*I[8]*I[25] + I[8]*I[26] - 1/2*I[10]*I[21] + 1/2*I[12]*I[13]);

    _<x,y,z> := PolynomialRing(Parent(I[1]), 3);
    Q := a11*x^2+2*a12*x*y+2*a13*x*z+a22*y^2+2*a23*y*z+a33*z^2;
    V := v11*x^2+2*v12*x*y+2*v13*x*z+v22*y^2+2*v23*y*z+v33*z^2;
    F := f111*x^3+3*f112*x^2*y+3*f113*x^2*z+3*f122*x*y^2+6*f123*x*y*z+3*f133*x*z^2+f222*y^3+3*f223*y^2*z+3*f233*y*z^2+f333*z^3;

    if minimize then
        Qphi, phi := MinimalModel(Conic(ProjectiveSpace(Parent(Q)), Q));
        Q := DefiningPolynomial(Qphi);
        F := Evaluate(F, DefiningPolynomials(phi));
        V := Evaluate(V, DefiningPolynomials(phi));
    end if;

    ct := GCD([Denominator(c) : c in Coefficients(V) cat Coefficients(F)]) /
            GCD([Numerator(c) : c in Coefficients(V) cat Coefficients(F)]);
    F *:= ct;
    V *:= ct;

    phi := ConicParametrization(Q);
    F := Evaluate(F, DefiningPolynomials(phi));
    V := Evaluate(V, DefiningPolynomials(phi));

    R := Parent(V);
    _<X,Y,Z,T> := PolynomialRing(BaseRing(R), 4);

    Coef_V := [MonomialCoefficient(V, (R.2)^i*(R.1)^(4-i)) : i in [0..4]];
    Coef_F := [MonomialCoefficient(F, (R.2)^i*(R.1)^(6-i)) : i in [0..6]];
    Q_final := X*Z-Y^2;
    E_final := ct*T^3+(Coef_V[1]*X^2+Coef_V[2]*X*Y+Coef_V[3]*X*Z+Coef_V[4]*Y*Z+Coef_V[5]*Z^2)*T+(Coef_F[1]*X^3+Coef_F[2]*X^2*Y+Coef_F[3]*X^2*Z+Coef_F[4]*X*Y*Z+Coef_F[5]*X*Z^2+Coef_F[6]*Y*Z^2+Coef_F[7]*Z^3);
    
    return Q_final, E_final;
end function;

intrinsic ReconstructionGenus4(I::SeqEnum : minimize := true) -> RngMPolElt, RngMPolElt
{Given a generic tuple of invariants for non-hyperelliptic curves of genus 4, returns a degree 2 polynomial Q and a degree 3 polynomial E in 4 variables,
such that V(Q,E) is a non-hyp genus 4 curve with said invariants (up to scaling). 
Both cases rank 3 and 4 are covered. 
If minimize is set to true, the returned forms Q, E are transformed as to have smaller coefficients.
}
	if #I eq 60 then
		Q, E := ReconstructionGenus4Rank3(I : minimize := minimize);
		return Q, E;
	elif #I eq 65 then 
		Q, E := ReconstructionGenus4Rank4(I : minimize := minimize);
		return Q, E;
    end if;
end intrinsic;


