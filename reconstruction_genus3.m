function Transvectant(f, g, h, d)
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
                        Tfgh +:= (-1)^(l+m+n)*Factorial(d)/(Factorial(i)*Factorial(j)*Factorial(k)*Factorial(l)*Factorial(m)*Factorial(n))*(Sf[i+l+1][j+m+1][k+n+1]*Sg[k+m+1][i+n+1][j+l+1]*Sh[j+n+1][k+l+1][i+m+1]);
                    end for;
                end for;
            end for;
		end for;
    end for;

	return Tfgh;
end function;

/* juste un exemple
S<a,b,c,f,g,h> := PolynomialRing(Rationals(), 6);
_<x,y,z> := PolynomialRing(S, [1,1,1]);

f0 := a*x^2+b*y^2+c*z^2+2*f*y*z+2*g*x*z+2*h*x*y;
Transvectant(f0, f0, f0, 2);
*/


/*
X1 := 3*x+5*y-7*z;
X2 := 4*x+12*y-z;
X3 := -x+6*y+2*z;

X := [X1, X2, X3];

Y1 := X1;
Y2 := X2;
Y3 := X3;

Y := [Y1, Y2, Y3];
f0 := Random(-5,5)*x^4+Random(-5,5)*x^3*y+Random(-5,5)*x^3*z+Random(-5,5)*y^3*x+Random(-5,5)*y^3*z+Random(-5,5)*z^3*x+Random(-5,5)*z^3*y+Random(-5,5)*x^2*y^2+Random(-5,5)*x^2*z^2+Random(-5,5)*y^2*z^2+Random(-5,5)*x^2*y*z+Random(-5,5)*y^2*x*z+Random(-5,5)*z^2*x*y;

*/
R<a400, a310, a301, a220, a211, a202, a130, a121, a112, a103, a040, a031, a022, a013, a004, a, b, c, d, e, f, g, h, i> := PolynomialRing(Rationals(), 24);
S<X,Y,Z> := PolynomialRing(R, 3);

f0 := a400*X^4+a310*X^3*Y+a301*X^3*Z+a220*X^2*Y^2+a211*X^2*Y*Y+a202*X^2*Z^2+a130*X*Y^3+a121*X*Y^2*Z+a112*X*Y*Z^2+a103*X*Z^3+a040*Y^4+a031*Y^3*Z+a022*Y^2*Z^2+a013*Y*Z^3+a004*Z^3;
f1 := Evaluate(f0, [a*X+b*Y+c*Z, d*X+e*Y+f*Z, g*X+h*Y+i*Z]);
-1/20736*Transvectant(f0, f0, f0, 4);
Factorization(Transvectant(f1, f1, f1, 4));




function reconstruction(f0)
    C36 := Transvectant(f0, f0, f0, 2);   
    C52 := Transvectant(C36, f0, f0, 4);
    C112 := Transvectant(C52, C52, f0, 2);
    C145 := Transvectant(C36, f0, C52^2, 3);
    C251 := Transvectant(C145, f0, C52^2, 4);
    C311 := Transvectant(C145, f0, C52*C112, 4);
    C371 := Transvectant(C145, f0, C112^2, 4);

    if Transvectant(C251, C311, C371, 1) eq 0 then
        return "Problem, not basis";
    end if;

    X := [C251, C311, C371];
    Y := X;

    f1 := Transvectant(f0, X[2]^4, X[3]^4, 4)*Y[1]^4+
    Transvectant(f0, X[3]^4, X[1]^4, 4)*Y[2]^4+
    Transvectant(f0, X[1]^4, X[2]^4, 4)*Y[3]^4+
    16*Transvectant(f0, X[2]^3*X[3], X[3]^3*X[1], 4)*Y[1]^3*Y[2]+
    16*Transvectant(f0, X[2]^3*X[1], X[3]^3*X[2], 4)*Y[1]^3*Y[3]+
    16*Transvectant(f0, X[3]^3*X[2], X[1]^3*X[3], 4)*Y[2]^3*Y[1]+
    16*Transvectant(f0, X[3]^3*X[1], X[1]^3*X[2], 4)*Y[2]^3*Y[3]+
    16*Transvectant(f0, X[1]^3*X[2], X[2]^3*X[3], 4)*Y[3]^3*Y[1]+
    16*Transvectant(f0, X[1]^3*X[3], X[2]^3*X[1], 4)*Y[3]^3*Y[2]+
    36*Transvectant(f0, X[2]^2*X[3]^2, X[3]^2*X[1]^2, 4)*Y[1]^2*Y[2]^2+
    36*Transvectant(f0, X[2]^2*X[1]^2, X[3]^2*X[2]^2, 4)*Y[1]^2*Y[3]^2+
    36*Transvectant(f0, X[3]^2*X[1]^2, X[1]^2*X[2]^2, 4)*Y[2]^2*Y[3]^2+
    -144*Transvectant(f0, X[2]^2*X[3]*X[1], X[3]^2*X[1]*X[2], 4)*Y[1]^2*Y[2]*Y[3]+
    -144*Transvectant(f0, X[3]^2*X[2]*X[1], X[1]^2*X[3]*X[2], 4)*Y[2]^2*Y[1]*Y[3]+
    -144*Transvectant(f0, X[1]^2*X[2]*X[3], X[2]^2*X[3]*X[1], 4)*Y[3]^2*Y[1]*Y[2];

    return f1/(Transvectant(X[1], X[2], X[3], 1)^4*2^9*3^3);
end function;

_<x,y,z> := PolynomialRing(GF(97), [1,1,1]);

c := 0;
for i in [1..1000] do
    f0 := Random(-5,5)*x^4+Random(-5,5)*x^3*y+Random(-5,5)*x^3*z+Random(-5,5)*y^3*x+Random(-5,5)*y^3*z+Random(-5,5)*z^3*x+Random(-5,5)*z^3*y+Random(-5,5)*x^2*y^2+Random(-5,5)*x^2*z^2+Random(-5,5)*y^2*z^2+Random(-5,5)*x^2*y*z+Random(-5,5)*y^2*x*z+Random(-5,5)*z^2*x*y;
    f1 := reconstruction(f0);
    if Type(f1) eq Type(f0) and f1 eq f0 then
        c +:= 1; 
    end if;
end for;
c;




///////////////////////////////////////////////////////////////////////////////
function Sum(L1, L2)
    res := 0;
    for i in [1..#L1] do
        res +:= L1[i]*L2[i];
    end for;
    return res;
end function;

S<b1, b2, b3> := PolynomialRing(Rationals(), 3);
R<x,y,z> := PolynomialRing(Rationals(), [1,1,1]);

function DualBasis(X)
    res := [];
    n := #X;
    M := [[[Evaluate(Transvectant(X[i], X[j], X[k], 1), [0,0,0]) : k in [1..n]] : j in [1..n]] : i in [1..n]];
    for i in [1..n] do
        b := [Random(-5,5) : i in [1..n]];
        N := Matrix([[Sum(M[j][k], b) : j in [1..n]] : k in [1..n]]);
        N;
        V := Vector(Rationals(), [0 : i in [1..n]]);
        V[i] := 1;
        V;
        a := Solution(N, V);
        Append(~res, [a, b]);
    end for;
    return res;
end function;

DualBasis([x, y, z]);


function Transvectant(f, g, h, d)
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
                        Tfgh +:= (-1)^(l+m+n)*Factorial(d)/(Factorial(i)*Factorial(j)*Factorial(k)*Factorial(l)*Factorial(m)*Factorial(n))*(Sf[i+l+1][j+m+1][k+n+1]*Sg[k+m+1][i+n+1][j+l+1]*Sh[j+n+1][k+l+1][i+m+1]);
                    end for;
                end for;
            end for;
		end for;
    end for;

	return Tfgh;
end function;


S<a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10> := PolynomialRing(Rationals(), 20);
R<x,y,z> := PolynomialRing(S, [1,1,1]);

function Sum(L1, L2)
    res := 0;
    for i in [1..#L1] do
        res +:= L1[i]*L2[i];
    end for;
    return res;
end function;

function DualBasis(X)
    res := [];
    n := #X;
    a := [a1,a2,a3,a4,a5,a6];
    //M := [Evaluate(Transvectant(X[i], Sum(X, a), Sum(X, b), 3), [0,0,0]): i in [1..n]];
    Matrix([[Transvectant(X[1], X[i], X[j], 2) : j in [1..n]] : i in [1..n]]);
    //M1 := Matrix([[Evaluate(Transvectant(X[i], X[j], Sum(X, b), 3), [0,0,0]) : j in [1..n]] : i in [1..n]]);
    L := Subsets(Set([1..n]));
    for l in L do 
        b := [0] cat [Random(-5,5) : i in [2..n]];
        for el in l do
            b[el] := 0;
        end for;
        M := [Evaluate(Transvectant(X[i], Sum(X, a), Sum(X, b), 3), [0,0,0]): i in [1..n]];
        I := Ideal([M[1]-1] cat M[2..n]);
        if GroebnerBasis(I) ne [1] then
           b;
           GroebnerBasis(I);
        end if;
    end for;
    //I2 := Ideal([M[1], M[2]-1, M[3]]);
    //I3 := Ideal([M[1], M[2], M[3]-1]);
    //I1 := Ideal([M[1]-1, M[2], M[3], M[4], M[5], M[6]]);
    //I2 := Ideal([M[1], M[2]-1, M[3], M[4], M[5], M[6]]);
    //I3 := Ideal([M[1], M[2], M[3]-1, M[4], M[5], M[6]]);

    return GroebnerBasis(I);
end function;

/*
P := Matrix(R, [[Random(-5,5) : i in [1..3]] : j in [1..3]]);
L := Transpose(Matrix([[x, y, z]]));
L := P*L;
L1 := [L[i][1] : i in [1..3]];
DualBasis(L1);
*/


R<a,b,c,d,e,g> := PolynomialRing(Rationals(), 6);
S<x,y,z> := PolynomialRing(R, 3);

X := [x,y,z];
P := Matrix(S, [[Random(-5,5) : i in [1..3]] : j in [1..3]]);
L := Transpose(Matrix([X]));
L := P*L;
L1 := [L[i][1] : i in [1..3]];
L2 := P^(-1)*L;
L2 := [L2[i][1] : i in [1..3]];

f := a*x+b*y+c*z;
X := L1;
Transvectant(Evaluate(f, L2), X[1], X[2], 1)*X[3]+
Transvectant(Evaluate(f, L2), X[2], X[3], 1)*X[1]+
Transvectant(Evaluate(f, L2), X[3], X[1], 1)*X[2];

f := a*x^2+b*x*y+c*x*z+d*y^2+e*y*z+g*z^2;

X := L1;
S1 := Transvectant(Evaluate(f, L2), X[2]^2, X[3]^2, 2)*X[1]^2;
S2 := Transvectant(Evaluate(f, L2), X[3]^2, X[1]^2, 2)*X[2]^2;
S3 := Transvectant(Evaluate(f, L2), X[1]^2, X[2]^2, 2)*X[3]^2;
S4 := 4*Transvectant(Evaluate(f, L2), X[2]*X[3], X[3]*X[1], 2)*X[1]*X[2];
S5 := 4*Transvectant(Evaluate(f, L2), X[2]*X[1], X[3]*X[2], 2)*X[1]*X[3];
S6 := 4*Transvectant(Evaluate(f, L2), X[3]*X[1], X[1]*X[2], 2)*X[2]*X[3];

S1+S2+S3+S4+S5+S6;

X := [x^2, x*y, x*z, y^2, y*z, z^2];
S := 0;
for i in [1..6] do
    for j in [1..6] do
        for k in [1..6] do 
            for l in [1..6] do
                for m in [1..6] do
                    for n in [1..6] do
                        S +:= 
                    end for;
                end for;
            end for;
        end for;
    end for;
end for;

for i in [1..3] do
    for j in [1..3] do
        for k in [1..3] do 
            S +:= P[i][3]*P[j][1]*P[k][2]*Transvectant(f, X[i], X[j], 1)*X[k];
        end for;
    end for;
end for;


for i in [1..3] do
    for j in [1..3] do
        for k in [1..3] do 
            S +:= P[i][2]*P[j][3]*P[k][1]*Transvectant(f, X[i], X[j], 1)*X[k];
        end for;
    end for;
end for;


//L1 := [x^3+x^2*y+4*y^3, x^2*y-x^3, x^2*z, x*y^2, x*y*z, x*z^2, 7*x*y^2+y^3, y^2*z, y*z^2, z^3];
DualBasis(L1);

function Sum(L1, L2, k)
    res := 0;
    for i in [k..#L1] do
        res +:= L1[i]*L2[i];
    end for;
    return res;
end function;


function Sum2(L1, L2, k)
    res := 0;
    for i in [1..k] do
        res +:= L1[i]*L2[i];
    end for;
    return res;
end function;

function NewIndex(L, i0, j0)
    res := [];
    for l in L do
        if l eq i0 then
            Append(~res, j0);
        elif l eq j0 then
            Append(~res, i0);
        else
            Append(~res, l);
        end if;
    end for;
    return res;
end function;


function Swap(A, i0, j0)
    res := A;
    d := #Random(Keys(A));
    d;
    if IsPower(#Keys(A), d) then
        b, n := IsPower(#Keys(A), d);
    else 
        "n not defined";
    end if;
    list_indexes := [i : i in [0..n^d-1]];
    for l in list_indexes do
        index := Decomposition(l, n, d);
        new_index := NewIndex(index, i0, j0);
        new_index;
        res[index] := A[new_index];
    end for;
    return res;
end function;


function DualBasis(X)
    res := [];
    deg := Degree(X[1]);
    n := #X;
    M0 := [[[Evaluate(Transvectant(X[i], X[j], X[k], d), [0,0,0]) : k in [1..n]] : j in [1..n]] : i in [1..n]];
    for i in [1..n] do
        b := [Random(-1, 1) : j in [1..n]];
        try
            M := Swap();
            N := Transpose(Matrix([[Sum(M[j][k], b, k) : j in [1..n]] : k in [1..n]]));
            V := Vector(Rationals(), [0 : j in [1..n]]);
            V[i] := 1;
            a := Solution(N, V);
            a := [a[j] : j in [1..n]];
            Append(~res, [a, b]);
        catch e
            b := [Random(1, 20) : j in [1..n]];
            N := Transpose(Matrix([[Sum2(M[j][k], b, k) : j in [1..n]] : k in [1..n]]));
            V := Vector(Rationals(), [0 : j in [1..n]]);
            V[n] := 1;
            a := Solution(N, V);
            a := [a[j] : j in [1..n]];
            Append(~res, [a, b]);   
        end try;
    end for;
    return res;
end function;

R<x,y,z> := PolynomialRing(Rationals(), 3);

P := Matrix(R, [[Random(-5,5) : i in [1..10]] : j in [1..10]]);
L := Transpose(Matrix([[x^3, x^2*y, x^2*z, x*y^2, x*y*z, x*z^2, y^3, y^2*z, y*z^2, z^3]]));
L := P*L;
L1 := [L[i][1] : i in [1..10]];
L0 := DualBasis(L1);

for i in [1..#L1] do
    for j in [1..#L0] do
        i, j, Transvectant(L1[i], Sum(L0[j][2], L1, j), Sum2(L0[j][1], L1, j), 3);
    end for;
end for;

DualBasis([x^3, x^2*y, x^2*z, x*y^2, x*y*z, x*z^2, y^3, y^2*z, y*z^2, z^3]);



/*
function Decomposition(m, n, d)
    res := [];
    while #res lt d do
        Append(~res, (m mod n) +1);
        m := m div n;
    end while;
    return res;
end function;

function ConstructionTransvectantMatrix(X, d);
    n := #X;
    res := AssociativeArray();
    list_indexes := [i : i in [0..n^d-1]];
    for l in list_indexes do
        res[Decomposition(l, n, d)] := Decomposition(l, n, d);
    end for;
    return res;
end function;

function Evaluation(M, L)
    res_int := M;
    for l in L do
        res_int := res_int[l];
    end for;
    return res_int;
end function;

function Swap(A, i0, j0)
    res := A;
    d := #Random(Keys(A));
    d;
    if IsPower(#Keys(A), d) then
        b, n := IsPower(#Keys(A), d);
    else 
        "n not defined";
    end if;
    list_indexes := [i : i in [0..n^d-1]];
    for l in list_indexes do
        index := Decomposition(l, n, d);
        new_index := NewIndex(index, i0, j0);
        new_index;
        res[index] := A[new_index];
    end for;
    return res;
end function;

A := ConstructionTransvectantMatrix([1,2,3], 2);
A;
B := Swap(A, 1, 2);
for l in Keys(B) do
    l, A[l], B[l];"";
end for;
*/