import "Decomposition_genus_3/i11.m" : decompo11, L11, ind11;
import "Decomposition_genus_3/i12.m" : decompo12, L12, ind12;
import "Decomposition_genus_3/i13.m" : decompo13, L13, ind13;
import "Decomposition_genus_3/i21.m" : decompo21, L21, ind21;
import "Decomposition_genus_3/i22.m" : decompo22, L22, ind22;
import "Decomposition_genus_3/i23.m" : decompo23, L23, ind23;
import "Decomposition_genus_3/i31.m" : decompo31, L31, ind31;
import "Decomposition_genus_3/i32.m" : decompo32, L32, ind32;
import "Decomposition_genus_3/i33.m" : decompo33, L33, ind33;

import "Decomposition_genus_3/i1111.m" : decompo1111, L1111, ind1111;
import "Decomposition_genus_3/i1112.m" : decompo1112, L1112, ind1112;
import "Decomposition_genus_3/i1113.m" : decompo1113, L1113, ind1113;
import "Decomposition_genus_3/i1122.m" : decompo1122, L1122, ind1122;
import "Decomposition_genus_3/i1123.m" : decompo1123, L1123, ind1123;
import "Decomposition_genus_3/i1133.m" : decompo1133, L1133, ind1133;
import "Decomposition_genus_3/i1222.m" : decompo1222, L1222, ind1222;
import "Decomposition_genus_3/i1223.m" : decompo1223, L1223, ind1223;
import "Decomposition_genus_3/i1233.m" : decompo1233, L1233, ind1233;
import "Decomposition_genus_3/i1333.m" : decompo1333, L1333, ind1333;
import "Decomposition_genus_3/i2222.m" : decompo2222, L2222, ind2222;
import "Decomposition_genus_3/i2223.m" : decompo2223, L2223, ind2223;
import "Decomposition_genus_3/i2233.m" : decompo2233, L2233, ind2233;
import "Decomposition_genus_3/i2333.m" : decompo2333, L2333, ind2333;
import "Decomposition_genus_3/i3333.m" : decompo3333, L3333, ind3333;

function Prod(L1, L2)
    return &*[L2[i] : i in L1];
end function;

function evaluation(decompo, L, ind, inv)
    m := #L;
    n := #decompo;
    res := &+[decompo[i]*Prod(L[ind[i]], inv) : i in [1..n] | ind[i] le m];
    return res;
end function;

function ReconstructionFromDixmierOhnoBis(inv : Minimize := false)
    _<X,Y,Z> := PolynomialRing(Rationals(), 3);
    
    f := evaluation(decompo1111, L1111, ind1111, inv)*X^4+
    4*evaluation(decompo1112, L1112, ind1112, inv)*X^3*Y+
    4*evaluation(decompo1113, L1113, ind1113, inv)*X^3*Z+
    6*evaluation(decompo1122, L1122, ind1122, inv)*X^2*Y^2+
    12*evaluation(decompo1123, L1123, ind1123, inv)*X^2*Y*Z+
    6*evaluation(decompo1133, L1133, ind1133, inv)*X^2*Z^2+
    4*evaluation(decompo1222, L1222, ind1222, inv)*X*Y^3+
    12*evaluation(decompo1223, L1223, ind1223, inv)*X*Y^2*Z+
    12*evaluation(decompo1233, L1233, ind1233, inv)*X*Y*Z^2+
    4*evaluation(decompo1333, L1333, ind1333, inv)*X*Z^3+
    evaluation(decompo2222, L2222, ind2222, inv)*Y^4+
    4*evaluation(decompo2223, L2223, ind2223, inv)*Y^3*Z+
    6*evaluation(decompo2233, L2233, ind2233, inv)*Y^2*Z^2+
    4*evaluation(decompo2333, L2333, ind2333, inv)*Y*Z^3+
    evaluation(decompo3333, L3333, ind3333, inv)*Z^4;
    
    if Minimize then
        return MinRedTernaryForm(f);
    end if;
    
    return f;
end function;


// Random example

_<x,y,z> := PolynomialRing(Rationals(), 3);
f := Random(-5, 5)*x^4+Random(-5,5)*x^3*y+Random(-5,5)*x^3*z+Random(-5,5)*y^3*x+Random(-5,5)*y^3*z+Random(-5,5)*z^3*x+Random(-5,5)*z^3*y+Random(-5,5)*x^2*y^2+Random(-5,5)*x^2*z^2+Random(-5,5)*y^2*z^2+Random(-5,5)*x^2*y*z+Random(-5,5)*y^2*x*z+Random(-5,5)*z^2*x*y+Random(-5,5)*y^4+Random(-5,5)*z^4;

inv := DixmierOhnoInvariants(f);

f1 := ReconstructionFromDixmierOhnoBis(inv);
