////////////////////////////////////// reconstruction of genus 4 curves with biforms of bidegree (1,1)

function MatrixTransvectant(X)
    return Matrix(BaseRing(Parent(X[1])), [[Transvectant(X[i], X[j], 1, 1 : invariant := true) : i in [1..4]] : j in [1..4]]);
end function;

function IsBasis(X)
    M := MatrixTransvectant(X);
    det := Determinant(M);
    if det eq 0 then
        return false;
    else 
        return true;
    end if;
end function;

function DualBasis(X)
    res := [];
    M := MatrixTransvectant(X);
    N := Determinant(M)*M^(-1);
    for i in [1..4] do
        V := Vector(Rationals(), [0,0,0,0]);
        V[i] := 1;
        sol := V*Transpose(N);
        res_int := 0;
        for i in [1..4] do
            res_int +:= sol[i]*X[i];
        end for;
        Append(~res, res_int);
    end for;
    return res;
end function;

// Checks that the quadric vanishes
procedure Verif(X, X2)
    Matrix([[Transvectant(X[i], X2[j], 1, 1 : invariant := true) : j in [1..4]] : i in [1..4]]);

    S := 0;
    for i in [1..4] do
        S +:= X[i]*X2[i];
    end for;
    S;

    S := 0;
    for i in [1..4] do
        for j in [1..4] do
            S +:= Transvectant(X[i], X[j], 1, 1)*X2[i]*X2[j];
        end for;
    end for;
    S;
end procedure;

// given a bicubic form f, and a basis of the space of biforms of bidegree (1,1), returns the quadric and cubic reconstructed from that basis.
function QuadricCubic(f, X)
    RP3<W,Y,Z,T> := PolynomialRing(Rationals(), 4);
    Q := 0;
    C := 0;
    for i in [1..4] do
        for j in [1..4] do
            for k in [1..4] do
                C +:= Transvectant(f, X[i]*X[j]*X[k], 3, 3 : invariant := true)*RP3.i*RP3.j*RP3.k;
            end for;
        Q +:= Transvectant(X[i], X[j], 1, 1 : invariant := true)*RP3.i*RP3.j;
        end for;
    end for;
    return Q, C;
end function;


// Given invariants, reconstructs a quadric and a cubic given those invariants. Only works when the covariants c3, c51, c52, c53 form a basis of the space of biforms of bidegree (1,1). We can use other covariants to work in more generality
function ReconstructionFromInvariants(L)
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
	c3c71 := -729/116*L[1]^2*L[5]+183239/268656*L[1]^2*L[4]+1/4*L[2]*L[4]-1407397/1432832*L[1]*L[2]*L[3]+112707525/5731328*L[1]^5-61121187/2865664*L[1]^3*L[3]+1407397/1432832*L[1]^3*L[2]+5645/1544*L[1]*L[6]+14/3*L[8]
			+729/116*L[3]*L[5]+9534849/5731328*L[1]*L[3]^2-99/464*L[3]*L[4]+60993/772*L[1]*L[12]-215/1158*L[1]*L[13]+2673/928*L[3]*L[11]+5/8*L[2]*L[11]-715391785/34567072*L[1]^2*L[11]-6378996231/19071488*L[14]+11016429/595984*L[15]
			+6312303/9535744*L[16]-25465/74498*L[17]-13941/595984*L[18]-8343/19071488*L[19];
	c51c71 := L[20];
	c52c71 := L[21];
	c53c71 := L[22];
	c71c71 := L[29];

	M_3_51_52_53 := Matrix(Parent(L[1]), [[c3c3, c3c51, c3c52, c3c53], [c3c51, c51c51, c51c52, c51c53], [c3c52, c51c52, c52c52, c52c53], [c3c53, c51c53, c52c53, c53c53]]);
	M_3_51_52_71 := Matrix(Parent(L[1]), [[c3c3, c3c51, c3c52, c3c71], [c3c51, c51c51, c51c52, c51c71], [c3c52, c51c52, c52c52, c52c71], [c3c71, c51c71, c52c71, c71c71]]);
	M_3_51_53_71 := Matrix(Parent(L[1]), [[c3c3, c3c51, c3c53, c3c71], [c3c51, c51c51, c51c53, c51c71], [c3c53, c51c53, c53c53, c53c71], [c3c71, c51c71, c53c71, c71c71]]);
	M_3_52_53_71 := Matrix(Parent(L[1]), [[c3c3, c3c52, c3c53, c3c71], [c3c52, c52c52, c52c53, c52c71], [c3c53, c52c53, c53c53, c53c71], [c3c71, c52c71, c53c71, c71c71]]);
	M_51_52_53_71 := Matrix(Parent(L[1]), [[c51c51, c51c52, c51c53, c51c71], [c51c52, c52c52, c52c53, c52c71], [c51c53, c52c53, c53c53, c53c71], [c51c71, c52c71, c53c71, c71c71]]);

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
		C := 1/Determinant(M_3_51_52_53)^3*(fc3c3c3*X^3+3*fc3c3c51*X^2*Y+3*fc3c3c52*X^2*Z+3*fc3c3c53*X^2*T+3*fc3c51c51*X*Y^2+6*fc3c51c52*X*Y*Z+6*fc3c51c53*X*Y*T+3*fc3c52c52*X*Z^2+6*fc3c52c53*X*Z*T+3*fc3c53c53*X*T^2+fc51c51c51*Y^3+3*fc51c51c52*Y^2*Z+3*fc51c51c53*Y^2*T+3*fc51c52c52*Y*Z^2+6*fc51c52c53*Y*Z*T+3*fc51c53c53*Y*T^2+fc52c52c52*Z^3+3*fc52c52c53*Z^2*T+3*fc52c53c53*Z*T^2+fc53c53c53*T^3); 
	
	elif Determinant(M_3_51_52_71) ne 0 then
		fc3c3c3 := L[8];
		fc3c3c51 := L[26];
		fc3c3c52 := L[27];
		fc3c3c71 := L[40];
		fc3c51c51 := L[34];
		fc3c51c52 := L[35];
		fc3c51c71 := L[41];
		fc3c52c52 := L[37];
		fc3c52c71 := L[42];
		fc3c71c71 := L[55];
		fc51c51c51 := L[45];
		fc51c51c52 := L[46];
		fc51c51c71 := L[56];
		fc51c52c52 := L[48];
		fc51c52c71 := L[57];
		fc51c71c71 := 0;
		fc52c52c52 := L[51];
		fc52c52c71 := L[59];
		fc52c71c71 := 0;
		fc71c71c71 := 0;

		Q := c3c3*X^2+2*c3c51*X*Y+2*c3c52*X*Z+2*c3c71*X*T+c51c51*Y^2+2*c51c52*Y*Z+2*c51c71*Y*T+c52c52*Z^2+2*c52c71*Z*T+c71c71*T^2;
		C := fc3c3c3*X^3+3*fc3c3c51*X^2*Y+3*fc3c3c52*X^2*Z+3*fc3c3c71*X^2*T+3*fc3c51c51*X*Y^2+6*fc3c51c52*X*Y*Z+6*fc3c51c71*X*Y*T+3*fc3c52c52*X*Z^2+6*fc3c52c71*X*Z*T+3*fc3c71c71*X*T^2+fc51c51c51*Y^3+3*fc51c51c52*Y^2*Z+3*fc51c51c71*Y^2*T+3*fc51c52c52*Y*Z^2+6*fc51c52c71*Y*Z*T+3*fc51c71c71*Y*T^2+fc52c52c52*Z^3+3*fc52c52c71*Z^2*T+3*fc52c71c71*Z*T^2+fc71c71c71*T^3; 
	
	elif Determinant(M_3_51_53_71) ne 0 then
		fc3c3c3 := L[8];
		fc3c3c51 := L[26];
		fc3c3c53 := L[28];
		fc3c3c71 := L[40];
		fc3c51c51 := L[34];
		fc3c51c53 := L[36];
		fc3c51c71 := L[41];
		fc3c53c53 := L[39];
		fc3c53c71 := L[43];
		fc3c71c71 := L[55];
		fc51c51c51 := L[45];
		fc51c51c53 := L[47];
		fc51c51c71 := L[56];
		fc51c53c53 := L[50];
		fc51c53c71 := L[58];
		fc51c71c71 := 0;
		fc53c53c53 := L[54];
		fc53c53c71 := L[61];
		fc53c71c71 := 0;
		fc71c71c71 := 0;

		Q := c3c3*X^2+2*c3c51*X*Y+2*c3c53*X*Z+2*c3c71*X*T+c51c51*Y^2+2*c51c53*Y*Z+2*c51c71*Y*T+c53c53*Z^2+2*c53c71*Z*T+c71c71*T^2;
		C := fc3c3c3*X^3+3*fc3c3c51*X^2*Y+3*fc3c3c53*X^2*Z+3*fc3c3c71*X^2*T+3*fc3c51c51*X*Y^2+6*fc3c51c53*X*Y*Z+6*fc3c51c71*X*Y*T+3*fc3c53c53*X*Z^2+6*fc3c53c71*X*Z*T+3*fc3c71c71*X*T^2+fc51c51c51*Y^3+3*fc51c51c53*Y^2*Z+3*fc51c51c71*Y^2*T+3*fc51c53c53*Y*Z^2+6*fc51c53c71*Y*Z*T+3*fc51c71c71*Y*T^2+fc53c53c53*Z^3+3*fc53c53c71*Z^2*T+3*fc53c71c71*Z*T^2+fc71c71c71*T^3; 
	
	elif Determinant(M_3_52_53_71) ne 0 then
		fc3c3c3 := L[8];
		fc3c3c52 := L[27];
		fc3c3c53 := L[28];
		fc3c3c71 := L[40];
		fc3c52c52 := L[37];
		fc3c52c53 := L[38];
		fc3c52c71 := L[42];
		fc3c53c53 := L[39];
		fc3c53c71 := L[43];
		fc3c71c71 := L[55];
		fc52c52c52 := L[51];
		fc52c52c53 := L[52];
		fc52c52c71 := L[59];
		fc52c53c53 := L[53];
		fc52c53c71 := L[60];
		fc52c71c71 := 0;
		fc53c53c53 := L[54];
		fc53c53c71 := L[61];
		fc53c71c71 := 0;
		fc71c71c71 := 0;

		Q := c3c3*X^2+2*c3c52*X*Y+2*c3c53*X*Z+2*c3c71*X*T+c52c52*Y^2+2*c52c53*Y*Z+2*c52c71*Y*T+c53c53*Z^2+2*c53c71*Z*T+c71c71*T^2;
		C := (fc3c3c3*X^3+3*fc3c3c52*X^2*Y+3*fc3c3c53*X^2*Z+3*fc3c3c71*X^2*T+3*fc3c52c52*X*Y^2+6*fc3c52c53*X*Y*Z+6*fc3c52c71*X*Y*T+3*fc3c53c53*X*Z^2+6*fc3c53c71*X*Z*T+3*fc3c71c71*X*T^2+fc52c52c52*Y^3+3*fc52c52c53*Y^2*Z+3*fc52c52c71*Y^2*T+3*fc52c53c53*Y*Z^2+6*fc52c53c71*Y*Z*T+3*fc52c71c71*Y*T^2+fc53c53c53*Z^3+3*fc53c53c71*Z^2*T+3*fc53c71c71*Z*T^2+fc71c71c71*T^3); 
	
	elif Determinant(M_51_52_53_71) ne 0 then
		fc51c51c51 := L[45];
		fc51c51c52 := L[46];
		fc51c51c53 := L[47];
		fc51c51c71 := L[56];
		fc51c52c52 := L[48];
		fc51c52c53 := L[49];
		fc51c52c71 := L[57];
		fc51c53c53 := L[50];
		fc51c53c71 := L[58];
		fc51c71c71 := 0;
		fc52c52c52 := L[51];
		fc52c52c53 := L[52];
		fc52c52c71 := L[59];
		fc52c53c53 := L[53];
		fc52c53c71 := L[60];
		fc52c71c71 := 0;
		fc53c53c53 := L[54];
		fc53c53c71 := L[61];
		fc53c71c71 := 0;
		fc71c71c71 := 0;

		Q := c51c51*X^2+2*c51c52*X*Y+2*c51c53*X*Z+2*c51c71*X*T+c52c52*Y^2+2*c52c53*Y*Z+2*c52c71*Y*T+c53c53*Z^2+2*c53c71*Z*T+c71c71*T^2;
		C := fc51c51c51*X^3+3*fc51c51c52*X^2*Y+3*fc51c51c53*X^2*Z+3*fc51c51c71*X^2*T+3*fc51c52c52*X*Y^2+6*fc51c52c53*X*Y*Z+6*fc51c52c71*X*Y*T+3*fc51c53c53*X*Z^2+6*fc51c53c71*X*Z*T+3*fc51c71c71*X*T^2+fc52c52c52*Y^3+3*fc52c52c53*Y^2*Z+3*fc52c52c71*Y^2*T+3*fc52c53c53*Y*Z^2+6*fc52c53c71*Y*Z*T+3*fc52c71c71*Y*T^2+fc53c53c53*Z^3+3*fc53c53c71*Z^2*T+3*fc53c71c71*Z*T^2+fc71c71c71*T^3; 
	else
		Q := 0;
		C := 0;
	end if;
	
	return Q, C;
end function;



function InvariantsGenus4CurvesRank4(f : normalize := false)
	K := BaseRing(Parent(f));
	
    GCD_hsop := [288, 12288, 746496, 12582912, 1741425868800, 19327352832, 764411904, 144, 570630428688384, 4076863488];
	GCD_others := [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 144, 144, 144, 1, 1, 1, 1, 1, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144 ];

	Jac := Transvectant(f, f, 1, 1);
	H := Transvectant(f, f, 2, 2);

	// Covariants
	// Degree 3
	c31 := ExactQuotient(Transvectant(H, f, 2, 2), 1536);
	c331 := Transvectant(Jac, f, 2, 2);
	c332 := Transvectant(H, f, 1, 1);

	// Degree 4
	c421 := Transvectant(H, H, 1, 1);
	c422 := Transvectant(c31, f, 1, 1);
	c423 := Transvectant(c332, f, 2, 2);

	c441 := Transvectant(c332, f, 1, 1);
	c442 := Transvectant(Transvectant(Jac, f, 1, 1), f, 2, 2);

	// Degree 5
	c511 := ExactQuotient(Transvectant(c422, f, 2, 2), 144);
	c512 := ExactQuotient(Transvectant(c441, f, 3, 3), 995328);
	c513 := ExactQuotient(Transvectant(c442, f, 3, 3), 1492992);
	
	c531 := Transvectant(c422, f, 1, 1);
	c532 := Transvectant(c423, f, 1, 1);
	c533 := Transvectant(Transvectant(f^3, f, 3, 3), f, 3, 3);

	// Degree 6
	c621 := Transvectant(c531, f, 2, 2);
	c622 := Transvectant(c532, f, 2, 2);
	c623 := Transvectant(c511, f, 1, 1);

	// Degree 7
	c711 := ExactQuotient(Transvectant(c621, f, 2, 2), 9216);
	c712 := ExactQuotient(Transvectant(c511, Transvectant(f, f, 2, 2), 1, 1), 32);
	c713 := ExactQuotient(Transvectant(c512, Transvectant(f, f, 2, 2), 1, 1), 96);
    
	c731 := Transvectant(c622, f, 1, 1);
	c732 := Transvectant(c623, f, 1, 1);

	// Degree 8 
	c821 := Transvectant(c711, f, 1, 1);
	c822 := Transvectant(c732, f, 2, 2);

	c84 := Transvectant(c731, f, 1, 1);

	// Degree 9
	c91 := Transvectant(c822, f, 2, 2);
	
	c931 := Transvectant(c821, f, 1, 1);
	c932 := Transvectant(c84, f, 2, 2);

	// Degree 10
	c102 := Transvectant(c91, f, 1, 1);

	// Degree 11
	c111 := Transvectant(c102, f, 2, 2);

	c113 := Transvectant(c932*f, f, 3, 3);

	// Invariants
	// HSOP
	I2 := Transvectant(f, f, 3, 3 : invariant := true);//288
	I41 := Transvectant(H, H, 2, 2 : invariant := true);//12288
	I42 := Transvectant(c331, f, 3, 3 : invariant := true);//746496
	I61 := Transvectant(H, c421, 2, 2 : invariant := true);//12582912
	I62 := Transvectant(c533, f, 3, 3 : invariant :=  true);//1741425868800
	I81 := Transvectant(c421, c421, 2, 2 : invariant := true);//19327352832
	I82 := Transvectant(c731, f, 3, 3 : invariant := true);//764411904
	I10 := Transvectant(f, c31^3, 3, 3 : invariant := true);//432
	I12 := Transvectant(c113, f, 3, 3 : invariant := true);//570630428688384
	I14 := Transvectant(c111*H, f, 3, 3 : invariant := true);//4076863488
	invHSOP := [K | I2,I41,I42,I61,I62,I81,I82,I10,I12,I14];

	// Degree 6
	j61 := Transvectant(c31, c31, 1, 1 : invariant := true);//2
	inv6 := [K | j61];

	// Degree 8
	j81 := Transvectant(c31, c511, 1, 1 : invariant := true);//2
	j82 := Transvectant(c31, c512, 1, 1 : invariant := true);//2
	inv8 := [K | j81,j82];
	
	// Degree 10
	j101 := Transvectant(c511, c511, 1, 1 : invariant := true);//2
    j102 := Transvectant(c511, c512, 1, 1 : invariant := true);//2
    j103 := Transvectant(c511, c513, 1, 1 : invariant := true);//2
    j104 := Transvectant(c512, c512, 1, 1 : invariant := true);//2
    j105 := Transvectant(c512, c513, 1, 1 : invariant := true);//2
    j106 := Transvectant(c513, c513, 1, 1 : invariant := true);//2
	inv10 := [K | j101,j102,j103,j104,j105,j106];

	// Degree 12	
    j121 := Transvectant(c711, c511, 1, 1 : invariant := true);//1
    j122 := Transvectant(c711, c512, 1, 1 : invariant := true);//2
    j123 := Transvectant(c711, c513, 1, 1 : invariant := true);//1
    j124 := Transvectant(c712, c511, 1, 1 : invariant := true);//2
    j125 := Transvectant(c712, c512, 1, 1 : invariant := true);//6
    j126 := Transvectant(c712, c513, 1, 1 : invariant := true);//2
    j127 := Transvectant(f, c511*c31^2, 3, 3 : invariant := true);//432
    j128 := Transvectant(f, c512*c31^2, 3, 3 : invariant := true);//432
    j129 := Transvectant(f, c513*c31^2, 3, 3 : invariant := true);//432
	inv12 := [K | j121,j122,j123,j124,j125,j126,j127,j128,j129];

	// Degree 14
	j141 := Transvectant(c711, c711, 1, 1 : invariant := true);//2
    j142 := Transvectant(c711, c712, 1, 1 : invariant := true);//1
	j143 := Transvectant(c711, c713, 1, 1 : invariant := true);//2
    j144 := Transvectant(c712, c713, 1, 1 : invariant := true);//2
	j145 := Transvectant(c713, c713, 1, 1 : invariant := true);//2
    j146 := Transvectant(f, c511*c511*c31, 3, 3 : invariant := true);//144
    j147 := Transvectant(f, c511*c512*c31, 3, 3 : invariant := true);//432
    j148 := Transvectant(f, c511*c513*c31, 3, 3 : invariant := true);//144
    j149 := Transvectant(f, c512*c512*c31, 3, 3 : invariant := true);//432
    j1410 := Transvectant(f, c512*c513*c31, 3, 3 : invariant := true);//432
    j1411 := Transvectant(f, c513*c513*c31, 3, 3 : invariant := true);//144
    j1412 := Transvectant(f, c711*c31*c31, 3, 3 : invariant := true);//432
	inv14 := [K | j141,j142,j143,j144,j145,j146,j147,j148,j149,j1410,j1411,j1412];

	// Degree 16
	j161 := Transvectant(f, c711*c511*c31, 3, 3 : invariant := true);//144
    j162 := Transvectant(f, c711*c512*c31, 3, 3 : invariant := true);//432
    j163 := Transvectant(f, c711*c513*c31, 3, 3 : invariant := true);//144
    j164 := Transvectant(f, c712*c511*c31, 3, 3 : invariant := true);//144
    j165 := Transvectant(f, c511*c511*c511, 3, 3 : invariant := true);//432
    j166 := Transvectant(f, c511*c511*c512, 3, 3 : invariant := true);//144
    j167 := Transvectant(f, c511*c511*c513, 3, 3 : invariant := true);//144
    j168 := Transvectant(f, c511*c512*c512, 3, 3 : invariant := true);//432
    j169 := Transvectant(f, c511*c512*c513, 3, 3 : invariant := true);//144
    j1610 := Transvectant(f, c511*c513*c513, 3, 3 : invariant := true);//144
    j1611 := Transvectant(f, c512*c512*c512, 3, 3 : invariant := true);//432
    j1612 := Transvectant(f, c512*c512*c513, 3, 3 : invariant := true);//432
    j1613 := Transvectant(f, c512*c513*c513, 3, 3 : invariant := true);//144
    j1614 := Transvectant(f, c513*c513*c513, 3, 3 : invariant := true);//432
    inv16 := [K | j161,j162,j163,j164,j165,j166,j167,j168,j169,j1610,j1611,j1612,j1613,j1614];

   
	// Degree 18
	j181 := Transvectant(f, c711*c711*c31, 3, 3 : invariant := true);//144
	j182 := Transvectant(f, c711*c511*c511, 3, 3 : invariant := true);//144
	j183 := Transvectant(f, c711*c511*c512, 3, 3 : invariant := true);//144
	j184 := Transvectant(f, c711*c511*c513, 3, 3 : invariant := true);//144
	j185 := Transvectant(f, c711*c512*c512, 3, 3 : invariant := true);//432
    j186 := Transvectant(f, c711*c512*c513, 3, 3 : invariant := true);//144
	j187 := Transvectant(f, c711*c513*c513, 3, 3 : invariant := true);//144
	j188 := Transvectant(f, c711*c712*c31, 3, 3 : invariant := true);//144
	j189 := Transvectant(f, c712*c712*c31, 3, 3 : invariant := true);//144
	j1810 := Transvectant(f, c712*c511*c512, 3, 3 : invariant := true);//144
	j1811 := Transvectant(f, c712*c512*c512, 3, 3 : invariant := true);//432
	inv18 := [K | j181,j182,j183,j184,j185,j186,j187,j188,j189,j1810,j1811];

	invOthers := inv6 cat inv8 cat inv10 cat inv12 cat inv14 cat inv16 cat inv18;

	if Type(K) eq RngInt then
		Inv := [ExactQuotient(invHSOP[i], GCD_hsop[i]) : i in [1..#invHSOP]] cat [ExactQuotient(invOthers[i], GCD_others[i]) : i in [1..#invOthers]];
	else
		Inv := [invHSOP[i]/GCD_hsop[i] : i in [1..#invHSOP]] cat [invOthers[i]/GCD_others[i] : i in [1..#invOthers]];
	end if;

	Wgt := [2,4,4,6,6,8,8,10,12,14,6,8,8,10,10,10,10,10,10,12,12,12,12,12,12,12,12,12,14,14,14,14,14,14,14,14,14,14,14,14,16,16,16,16,16,16,16,16,16,16,16,16,16,16,18,18,18,18,18,18,18,18,18,18,18];
	
	if normalize then
		return WPSNormalize(Wgt, Inv), Wgt;
	end if;

	return Inv, Wgt;
end function;



function Covariants(f)
	K := BaseRing(Parent(f));
	
	Jac := Transvectant(f, f, 1, 1);
	H := Transvectant(f, f, 2, 2);

	// Covariants
	// Degree 3
	c31 := ExactQuotient(Transvectant(H, f, 2, 2), 1536);

	c331 := Transvectant(Jac, f, 2, 2);
	c332 := Transvectant(H, f, 1, 1);

	// Degree 4
	c421 := Transvectant(H, H, 1, 1);
	c422 := Transvectant(c31, f, 1, 1);
	c423 := Transvectant(c332, f, 2, 2);

	c441 := Transvectant(c332, f, 1, 1);
	c442 := Transvectant(Transvectant(Jac, f, 1, 1), f, 2, 2);

	// Degree 5
	c511 := ExactQuotient(Transvectant(c422, f, 2, 2), 144);
	c512 := ExactQuotient(Transvectant(c441, f, 3, 3), 995328);
	c513 := ExactQuotient(Transvectant(c442, f, 3, 3), 1492992);
	
	c531 := Transvectant(c422, f, 1, 1);
	c532 := Transvectant(c423, f, 1, 1);
	c533 := Transvectant(Transvectant(f^3, f, 3, 3), f, 3, 3);

	// Degree 6
	c621 := Transvectant(c531, f, 2, 2);
	c622 := Transvectant(c532, f, 2, 2);
	c623 := Transvectant(c511, f, 1, 1);

	// Degree 7
	c711 := ExactQuotient(Transvectant(c621, f, 2, 2), 9216);
	c712 := ExactQuotient(Transvectant(c511, Transvectant(f, f, 2, 2), 1, 1), 32);
	c713 := ExactQuotient(Transvectant(c512, Transvectant(f, f, 2, 2), 1, 1), 96);
    
	c731 := Transvectant(c622, f, 1, 1);
	c732 := Transvectant(c623, f, 1, 1);

	// Degree 8 
	c821 := Transvectant(c711, f, 1, 1);
	c822 := Transvectant(c732, f, 2, 2);

	c84 := Transvectant(c731, f, 1, 1);

	// Degree 9
	c91 := Transvectant(c822, f, 2, 2);
	
	c931 := Transvectant(c821, f, 1, 1);
	c932 := Transvectant(c84, f, 2, 2);

	// Degree 10
	c102 := Transvectant(c91, f, 1, 1);

	// Degree 11
	c111 := Transvectant(c102, f, 2, 2);

	c113 := Transvectant(c932*f, f, 3, 3);

	return [c31, c511, c512, c513, c711, c712, c713, c91, c111];
end function;


/*
res := [];
K := Rationals();
R<x,y,u,v> := PolynomialRing(K, [1,1,1,1]);
for j in [1..30] do
	f := 0;
	for i in [0..3] do
		for j in [0..3] do
			f +:= Random(-10,10)*x^i*y^(3-i)*u^j*v^(3-j);
		end for;
	end for;
	f := f+Evaluate(f, [y,x,u,v]);
	f;
	L, Wgt := InvariantsGenus4CurvesRank4(f);
	L;
	ReconstructionFromInvariants(L);
end for;
*/
K := Rationals();
K1 := GF(11);
S<x,y,u,v> := PolynomialRing(K, 4);
S1<X,Y,U,V> := PolynomialRing(K1, 4);
for k in [1..10000] do
	f := 0;
	for i in [0..3] do
		for j in [0..3] do
			f +:= S!Random(-10,10)*x^i*y^(3-i)*u^j*v^(3-j);
		end for;
	end for;
	X := Covariants(S1!f);
	bool := false;
	for i1 in [1..9] do
		for i2 in [i1+1..9] do
			for i3 in [i2+1..9] do
				for i4 in [i3+1..9] do
					if IsBasis([X[i1], X[i2], X[i3], X[i4]]) then	
						bool := true;
					end if;
				end for;
			end for;
		end for;
	end for;
	if not bool then
		f;
		X;
		L0 := InvariantsGenus4CurvesRank4(f);
		L := ChangeUniverse(L0, K1);
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
		c3c71 := -729/116*L[1]^2*L[5]+183239/268656*L[1]^2*L[4]+1/4*L[2]*L[4]-1407397/1432832*L[1]*L[2]*L[3]+112707525/5731328*L[1]^5-61121187/2865664*L[1]^3*L[3]+1407397/1432832*L[1]^3*L[2]+5645/1544*L[1]*L[6]+14/3*L[8]
				+729/116*L[3]*L[5]+9534849/5731328*L[1]*L[3]^2-99/464*L[3]*L[4]+60993/772*L[1]*L[12]-215/1158*L[1]*L[13]+2673/928*L[3]*L[11]+5/8*L[2]*L[11]-715391785/34567072*L[1]^2*L[11]-6378996231/19071488*L[14]+11016429/595984*L[15]
				+6312303/9535744*L[16]-25465/74498*L[17]-13941/595984*L[18]-8343/19071488*L[19];
		c51c71 := L[20];
		c52c71 := L[21];
		c53c71 := L[22];
		c71c71 := L[29];

		M_3_51_52_53 := Matrix(Parent(L[1]), [[c3c3, c3c51, c3c52, c3c53], [c3c51, c51c51, c51c52, c51c53], [c3c52, c51c52, c52c52, c52c53], [c3c53, c51c53, c52c53, c53c53]]);
		M_3_51_52_71 := Matrix(Parent(L[1]), [[c3c3, c3c51, c3c52, c3c71], [c3c51, c51c51, c51c52, c51c71], [c3c52, c51c52, c52c52, c52c71], [c3c71, c51c71, c52c71, c71c71]]);
		M_3_51_53_71 := Matrix(Parent(L[1]), [[c3c3, c3c51, c3c53, c3c71], [c3c51, c51c51, c51c53, c51c71], [c3c53, c51c53, c53c53, c53c71], [c3c71, c51c71, c53c71, c71c71]]);
		M_3_52_53_71 := Matrix(Parent(L[1]), [[c3c3, c3c52, c3c53, c3c71], [c3c52, c52c52, c52c53, c52c71], [c3c53, c52c53, c53c53, c53c71], [c3c71, c52c71, c53c71, c71c71]]);
		M_51_52_53_71 := Matrix(Parent(L[1]), [[c51c51, c51c52, c51c53, c51c71], [c51c52, c52c52, c52c53, c52c71], [c51c53, c52c53, c53c53, c53c71], [c51c71, c52c71, c53c71, c71c71]]);

		//M_3_51_52_53;"";
		//M_3_51_52_71;"";
		//M_3_51_53_71;"";
		//M_3_52_53_71;"";
		//M_51_52_53_71;"";
	end if;
end for;

K := Rationals();
_<x,y,u,v> := PolynomialRing(K, 4);
f := 0;
for i in [0..3] do
	for j in [0..3] do
		f +:= Random(-10,10)*x^i*y^(3-i)*u^j*v^(3-j);
	end for;
end for;

L, Wgt := InvariantsGenus4CurvesRank4(f : normalize := true);
L := ChangeUniverse(L, K);
Q, C := ReconstructionFromInvariants(L);
L1 := InvariantsGenus4Curves(Q, C : normalize := true);
L1 := ChangeUniverse(L1, K);
L1 eq L;

P, t := NewBasis(Q);
f0 := CubicNewBasis(Q,C);

R<x, y, u, v> := PolynomialRing(BaseRing(Parent(f0)), 4);
f_bic := Evaluate(f0, [x*u, y*u, x*v, y*v]);
		
Inv, Wgt := InvariantsGenus4CurvesRank4(f_bic : normalize := true);
		
Inv eq ChangeUniverse(L, Parent(Inv[1]));

//L1 := InvariantsGenus4Curves(Q, C); 
//L eq L1;


K := Rationals();
R<x,y,u,v> := PolynomialRing(K, [1,1,1,1]);

f := 0;
for i in [0..3] do
	for j in [0..3] do
		f +:= Random(-10,10)*x^i*y^(3-i)*u^j*v^(3-j);
	end for;
end for;



R<x> := PolynomialRing(Rationals(), 23);
L := [R.i : i in [1..22]] cat [0,0,0,0,0,0,R.23];

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
c3c71 := -729/116*L[1]^2*L[5]+183239/268656*L[1]^2*L[4]+1/4*L[2]*L[4]-1407397/1432832*L[1]*L[2]*L[3]+112707525/5731328*L[1]^5-61121187/2865664*L[1]^3*L[3]+1407397/1432832*L[1]^3*L[2]+5645/1544*L[1]*L[6]+14/3*L[8]
			+729/116*L[3]*L[5]+9534849/5731328*L[1]*L[3]^2-99/464*L[3]*L[4]+60993/772*L[1]*L[12]-215/1158*L[1]*L[13]+2673/928*L[3]*L[11]+5/8*L[2]*L[11]-715391785/34567072*L[1]^2*L[11]-6378996231/19071488*L[14]+11016429/595984*L[15]
			+6312303/9535744*L[16]-25465/74498*L[17]-13941/595984*L[18]-8343/19071488*L[19];
c51c71 := L[20];
c52c71 := L[21];
c53c71 := L[22];
c71c71 := L[29];

M_3_51_52_53 := Matrix(Parent(L[1]), [[c3c3, c3c51, c3c52, c3c53], [c3c51, c51c51, c51c52, c51c53], [c3c52, c51c52, c52c52, c52c53], [c3c53, c51c53, c52c53, c53c53]]);
M_3_51_52_71 := Matrix(Parent(L[1]), [[c3c3, c3c51, c3c52, c3c71], [c3c51, c51c51, c51c52, c51c71], [c3c52, c51c52, c52c52, c52c71], [c3c71, c51c71, c52c71, c71c71]]);
M_3_51_53_71 := Matrix(Parent(L[1]), [[c3c3, c3c51, c3c53, c3c71], [c3c51, c51c51, c51c53, c51c71], [c3c53, c51c53, c53c53, c53c71], [c3c71, c51c71, c53c71, c71c71]]);
M_3_52_53_71 := Matrix(Parent(L[1]), [[c3c3, c3c52, c3c53, c3c71], [c3c52, c52c52, c52c53, c52c71], [c3c53, c52c53, c53c53, c53c71], [c3c71, c52c71, c53c71, c71c71]]);
M_51_52_53_71 := Matrix(Parent(L[1]), [[c51c51, c51c52, c51c53, c3c71], [c3c52, c52c52, c52c53, c52c71], [c3c53, c52c53, c53c53, c53c71], [c3c71, c52c71, c53c71, c71c71]]);

I := Ideal([Determinant(M_3_51_52_53), Determinant(M_3_51_52_71), Determinant(M_3_51_53_71), Determinant(M_3_52_53_71)]);
NormalForm(Determinant(M_51_52_53_71), I);

///////////////////////////// reconstruction with biforms of bidegree (2,2) (obsolete)
/*
function MatrixTransvectant(X)
    return Matrix(BaseRing(Parent(X[1])), [[Transvectant(X[i], X[j], 2, 2) : i in [1..9]] : j in [1..9]]);
end function;

function IsBasis(X)
    M := MatrixTransvectant(X);
    det := Determinant(M);
    if det eq 0 then
        return false;
    else 
        return true;
    end if;
end function;

function DualBasis(X)
    res := [];
    M := MatrixTransvectant(X);
    N := Determinant(M)*M^(-1);
    for i in [1..9] do
        V := Vector(Rationals(), [0,0,0,0,0,0,0,0,0]);
        V[i] := 1;
        sol := V*Transpose(N);
        res_int := 0;
        for i in [1..9] do
            res_int +:= sol[i]*X[i];
        end for;
        Append(~res, res_int);
    end for;
    return res;
end function;

procedure Verif(X, X2)
    Matrix([[Transvectant(X[i],X2[j],2,2) : j in [1..9]] : i in [1..9]]);

    S := 0;
    for i in [1..9] do
        S +:= X[i]*X2[i];
    end for;
    S;

    S := 0;
    for i in [1..9] do
        for j in [1..9] do
            S +:= Transvectant(X[i], X[j], 2, 2)*X2[i]*X2[j];
        end for;
    end for;
    S;
end procedure;


R<x,y,u,v> := PolynomialRing(Rationals(), 4);
P := Matrix(R, [[Random(1,10) : i in [1..9]] : j in [1..9]]);
X0 := Vector(R, [x^2*u^2, x^2*u*v, x^2*v^2, x*y*u^2, x*y*u*v, x*y*v^2, y^2*u^2, y^2*u*v, y^2*v^2]);
X := X0*Transpose(P);

X2 := DualBasis(X);
Verif(X, X2);
*/



// Given invariants, reconstructs a quadric and a cubic given those invariants. Only works when the covariants c3, c51, c52, c53 form a basis of the space of biforms of bidegree (1,1). We can use other covariants to work in more generality
function ReconstructionFromInvariants(L)
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
		C := 1/Determinant(M_3_51_52_53)^3*(fc3c3c3*X^3+3*fc3c3c51*X^2*Y+3*fc3c3c52*X^2*Z+3*fc3c3c53*X^2*T+3*fc3c51c51*X*Y^2+6*fc3c51c52*X*Y*Z+6*fc3c51c53*X*Y*T+3*fc3c52c52*X*Z^2+6*fc3c52c53*X*Z*T+3*fc3c53c53*X*T^2+fc51c51c51*Y^3+3*fc51c51c52*Y^2*Z+3*fc51c51c53*Y^2*T+3*fc51c52c52*Y*Z^2+6*fc51c52c53*Y*Z*T+3*fc51c53c53*Y*T^2+fc52c52c52*Z^3+3*fc52c52c53*Z^2*T+3*fc52c53c53*Z*T^2+fc53c53c53*T^3); 

		return Q, C;
	end if;
	
	return "Not implemented yet";
end function;
	