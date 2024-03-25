function MatrixTransvectant(X, d)
    return Matrix(BaseRing(Parent(X[1])), [[Evaluate(Transvectant(X[i], X[j], d), [0,0]) : i in [1..#X]] : j in [1..#X]]);
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

function DualBasis(X, d)
    res := [];
    M := MatrixTransvectant(X, d);
    M;
    N := Determinant(M)*M^(-1);
    for i in [1..#X] do
        V := Vector(Rationals(), [0 : i in [1..#X]]);
        V[i] := 1;
        sol := V*Transpose(N);
        res_int := 0;
        for i in [1..#X] do
            res_int +:= sol[i]*X[i];
        end for;
        Append(~res, res_int);
    end for;
    return res;
end function;

procedure Verif(X, X2)
    Matrix([[Transvectant(X[i], X2[j], 1) : j in [1..2]] : i in [1..2]]);

    S := 0;
    for i in [1..2] do
        S +:= X[i]*X2[i];
    end for;
    S;

    S := 0;
    for i in [1..2] do
        for j in [1..2] do
            S +:= Transvectant(X[i], X[j], 1)*X2[i]*X2[j];
        end for;
    end for;
    S;
end procedure;

function Convert(n, m)
    res := [];
    while n gt 0 do
        Insert(~res, 1, (n mod 2) + 1);
        n := n div 2;
    end while;
    lgt := #res;
    res := [1 : i in [1..m-lgt]] cat res; 
    return res;
end function;

function Product(X, L)
    S := 1;
    for i in [1..#L] do
        S *:= X[L[i]];
    end for;
    return S;
end function;

function Form(f, X)
    RP<S,T> := PolynomialRing(Rationals(), 2);
    n := Degree(f);
    C := 0;
    for i in [0..2^n-1] do
        L := Convert(i, n);
        C +:= Evaluate(Transvectant(f, Product(X, L), n), [0,0])*Product([S,T], L);
    end for;
    return C;
end function;

function Invariants_5(f)
    K := BaseRing(Parent(f));
    C22 := Transvectant(f, f, 4);
    C33 := Transvectant(C22, f, 2);
    C51 := Transvectant(f, C22^2, 4);
    C62 := Transvectant(C33, C33, 2);
    C71 := Transvectant(C51, C22, 1);
    C111 := Transvectant(C51, C62, 1);
    I4 := K!Evaluate(Transvectant(C22, C22, 2), [0,0]);
    I8 := K!Evaluate(Transvectant(C62, C22, 2), [0,0]);
    I12 := K!Evaluate(Transvectant(C51, C71, 1), [0,0]);
    I16 := K!Evaluate(Transvectant(C111, C51, 1), [0,0]);
    I18 := K!Evaluate(Transvectant(C71, C111, 1), [0,0]);
    return [I4, I8, I12, I18];
end function;

function Reconstruction(L)
    _<x,y> := PolynomialRing(Parent(L[1]), 2);
    I4 := L[1];
    I8 := L[2];
    I12 := L[3];
    I18 := L[4];
    I16 := -1/3*I4*I8+1/6*I4*I12+1/2*I8^2;
    if I12 ne 0 then
        i26 := 2/3*I4^2*I18-I8*I18;
        i28 := -1/9*I4^5*I8+1/18*I4^4*I12+1/3*I4^3*I8^2+5/12*I4^2*I8*I12-1/4*I4*I8^3-5/12*I4*I12^2-3/4*I8^2*I12;
        i30 := -1/3*I4^3*I18+1/2*I4*I8*I18+I12*I18;
        i32 := 1/18*I4^6*I8-1/36*I4^5*I12-1/6*I4^4*I8^2-3/8*I4^3*I8*I12+1/8*I4^2*I8^3+7/24*I4^2*I12^2+5/8*I4*I8^2*I12+1/2*I8*I12^2;
        i34 := 1/6*I4^4*I18-1/4*I4^2*I8*I18-I4*I12*I18;
        i36 := -1/36*I4^7*I8+1/72*I4^6*I12+1/12*I4^5*I8^2+13/48*I4^4*I8*I12-1/16*I4^3*I8^3-3/16*I4^3*I12^2-7/16*I4^2*I8^2*I12-1/2*I4*I8*I12^2+I12^3;
        return (i26*x^5+5*i28*x^4*y+10*i30*x^3*y^2+10*i32*x^2*y^3+5*i34*x*y^4+i36*y^5);
    elif I16 ne 0 then
        i26 := 2/3*I4^2*I18-I8*I18;
        i32 := -1/18*I4^4*I8^2+1/9*I4^3*I8*I12+1/6*I4^2*I8^3-1/24*I4^2*I12^2-1/6*I4*I8^2*I12-1/8*I8^4+1/6*I8*I12^2;
        i38 := -1/6*I8*I12*I18;
        i44 := -1/54*I4^5*I8^3+1/36*I4^4*I8^2*I12+1/18*I4^3*I8^4-1/72*I4^3*I8*I12^2+1/36*I4^2*I8^3*I12+1/432*I4^2*I12^3-1/24*I4*I8^5-1/12*I4*I8^2*I12^2-5/48*I8^4*I12+1/36*I8*I12^3;
        i50 := -1/9*I4^2*I8^3*I18+1/9*I4*I8^2*I12*I18+1/6*I8^4*I18-1/36*I8*I12^2*I18;
        i56 := -1/54*I4^6*I8^4+1/27*I4^5*I8^3*I12+1/9*I4^4*I8^5-1/36*I4^4*I8^2*I12^2-11/54*I4^3*I8^4*I12+1/108*I4^3*I8*I12^3-5/24*I4^2*I8^6+5/36*I4^2*I8^3*I12^2-1/864*I4^2*I12^4+2/9*I4*I8^5*I12-1/24*I4*I8^2*I12^3+1/8*I8^7-17/288*I8^4*I12^2+1/216*I8*I12^4;
        return (i26*x^5+5*i32*x^4*y+10*i38*x^3*y^2+10*i44*x^2*y^3+5*i50*x*y^4+i56*y^5);
    elif I8 eq 0 and I4 ne 0 then
        return x^5+y^5;
    elif I8 eq 0 and I4 eq 0 then
        return x^5;
    else
        return x^4*y+x*y^4;
    end if;
end function;

K := Rationals();
R<x,y> := PolynomialRing(K, 2);

c := 0;
time for i in [1..100000] do
    f0 := Random(-10,10)*x^5+Random(-10,10)*x^4*y+Random(-10,10)*x^3*y^2+Random(-10,10)*x^2*y^3+Random(-10,10)*x*y^4+Random(-10,10)*y^5;
    L0 := Invariants_5(f0);
    f1 := Reconstruction(L0);
    if WPSNormalize([4,8,12,18], Invariants_5(f0)) eq WPSNormalize([4,8,12,18], Invariants_5(f1)) then
    //if Invariants_5(f1) eq Invariants_5(f0) then
        c +:= 1;
    else 
        "pbm";
        Invariants_5(f0);
        Invariants_5(f1);
    end if;
end for;




///////////////////////////////////:::
function MatrixTransvectant(X, d)
    return Matrix(BaseRing(Parent(X[1])), [[Evaluate(Transvectant(X[i], X[j], d), [0,0]) : i in [1..#X]] : j in [1..#X]]);
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

function DualBasis(X, d)
    res := [];
    M := MatrixTransvectant(X, d);
    M;
    N := Determinant(M)*M^(-1);
    for i in [1..#X] do
        V := Vector(Rationals(), [0 : i in [1..#X]]);
        V[i] := 1;
        sol := V*Transpose(N);
        res_int := 0;
        for i in [1..#X] do
            res_int +:= sol[i]*X[i];
        end for;
        Append(~res, res_int);
    end for;
    return res;
end function;

K := Rationals();
R<x,y> := PolynomialRing(K, 2);
L := Transpose(Matrix([[x^2, x*y, y^2]]));
P := Matrix(R, [[Random(-5,5) : i in [1..3]] : j in [1..3]]);
//P := Matrix([[1,0,0],[0,1,0],[0,0,1]]);
L := P*L;
L1 := [L[i][1] : i in [1..3]];
X := L1;
X2 := DualBasis(X, 2);

for i in [1..3] do
    for j in [1..3] do
        Matrix([[Transvectant(X2[i]*X2[j], X[k]*X[l], 4) : k in [1..3]] : l in [1..3]]);
        Rank(Matrix([[Transvectant(X2[i]*X2[j], X[k]*X[l], 4) : k in [1..3]] : l in [1..3]]));
    end for;
end for;


