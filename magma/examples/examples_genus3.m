R<X,Y,Z> := PolynomialRing(Rationals(), 3);

f := &+[Random(-5,5)*m : m in MonomialsOfDegree(R, 4)];

DO, wgt := DixmierOhnoInvariants(f);

f1 := R!ReconstructionGenus3(DO);

IsIsomorphicTernaryQuartics(f, f1);
