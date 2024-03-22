function DOBis(f)
    inv := DixmierOhnoInvariants(f);
    if Characteristic(Universe(inv)) in {19, 47, 277, 523} then inv[3] +:= inv[4]; end if;
    return inv;
end function;

load "decompo.m";

intrinsic ReconstructionFromDixmierOhnoBis(inv::SeqEnum : 
    minimize := false) -> RngMPolElt
{
Reconstruct a ternary quartic from a generic tuple of Dixmier-Ohno invariants. 
It is not implemented for characteristic 2, 3, 5, 7 or 17.

If minimize is set to true, the coefficients of the ternary quartic are minimized using Elsenhans and Stoll's method MinRedTernaryForm.
}

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
        if minimize then
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
