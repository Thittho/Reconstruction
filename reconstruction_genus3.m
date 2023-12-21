function DOBis(f)
    inv := DixmierOhnoInvariants(f);
    if Characteristic(Universe(inv)) in {19, 47, 277, 523} then inv[3] +:= inv[4]; end if;
    return inv;
end function;

load "decompo.m";

function ReconstructionFromDixmierOhnoBis(inv : Minimize := false)

    if Characteristic(Universe(inv)) in {2, 3, 5, 7, 17} then "Not implemented for char 2, 3, 5, 7 or 17 "; return 0; end if;

    M := Matrix([[P11(inv), P12(inv), P13(inv)],
                 [P21(inv), P22(inv), P23(inv)],
                 [P31(inv), P32(inv), P33(inv)]]);
    _<X,Y,Z> := PolynomialRing(FieldOfFractions(Universe(inv)), 3);
    
    f := P1111(inv)*X^4+
        4*P1112(inv)*X^3*Y+
        4*P1113(inv)*X^3*Z+
        6*P1122(inv)*X^2*Y^2+
        12*P1123(inv)*X^2*Y*Z+
        6*P1133(inv)*X^2*Z^2+
        4*P1222(inv)*X*Y^3+
        12*P1223(inv)*X*Y^2*Z+
        12*P1233(inv)*X*Y*Z^2+
        4*P1333(inv)*X*Z^3+
        P2222(inv)*Y^4+
        4*P2223(inv)*Y^3*Z+
        6*P2233(inv)*Y^2*Z^2+
        4*P2333(inv)*Y*Z^3+
        P3333(inv)*Z^4;
    
    if Determinant(M) ne 0 then
        if Minimize then
            f1 := MinRedTernaryForm(f);
            return f1;
        end if;
        return f;
    end if;
    
    inv1 := DOBis(f);
    if DixmierOhnoInvariantsEqual(inv, inv1) then
        return f;
    end if;
    
    "Not a basis, not implemented yet.";
    
    return 0;
end function;

procedure Reconstruction(F, n : a := 10, b := 10^2)
    R<x,y,z> := PolynomialRing(F, 3);
    L4 := MonomialsOfDegree(R, 4);
    if IsFinite(F) then
        forms := [&+[Random(F)*l : l in L4] : i in [1..10]];
    else
        forms := [&+[Random(a,b)*l : l in L4] : i in [1..10]];
    end if;

    for f in forms do
        inv := DOBis(f);

        time f1 := ReconstructionFromDixmierOhnoBis(inv); // Reconstruction
        I3,I6,I9,J9,I12,J12,I15,J15,I18,J18,I21,J21,I27 := Explode(inv);
        
        if f1 ne 0 then
            f1 := R!(f1);
            inv2 := DOBis(f1);
            DixmierOhnoInvariantsEqual(inv, inv2); // Check that the reconstructed form has the same invariants
        end if;
    end for;
end procedure;

// Reconstruction over a big field
K := GF(NextPrime(10^200)^10);
Reconstruction(K, 5);

// Over a small field
K := GF(41);
Reconstruction(K, 5);

// Over Q
K := Rationals();
Reconstruction(K, 5 : a := 10^5, b := 10^10);

// Example given in your article
ReconstructionFromDixmierOhnoBis([0,0,0,0,-7*19/(2^14*3^8*5^2),0,-11*19/(2^17*3^10*5^2),0,7*19^2/(2^20*3^11*5^3),19^2/(2^20*3^11*5^3),-19^2*31/(2^24*3^13*5^5),-17*19^2/(2^21*3^12*5^5),-19^2*6553/(2^39*3^6*5^5*11)]);
PlaneQuarticFromDixmierOhnoInvariants([0,0,0,0,-7*19/(2^14*3^8*5^2),0,-11*19/(2^17*3^10*5^2),0,7*19^2/(2^20*3^11*5^3),19^2/(2^20*3^11*5^3),-19^2*31/(2^24*3^13*5^5),-17*19^2/(2^21*3^12*5^5),-19^2*6553/(2^39*3^6*5^5*11)] : minimize := false);
