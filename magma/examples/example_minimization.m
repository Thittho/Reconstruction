R4<X,Y,Z,T> := PolynomialRing(Rationals(), 4);
mons2 := MonomialsOfDegree(R4, 2);
mons3 := MonomialsOfDegree(R4, 3);

Q := &+[Random(10)*m: m in mons2];
Gamma := &+[Random(10)*m: m in mons3];
Inv, Wgt := InvariantsGenus4Curves(Q, Gamma);
Inv := WPSMinimize(Wgt, Inv);
Qbis, Gammabis := ReconstructionGenus4(Inv);
	
time Q1, Gamma1 := MinimizeG4(Qbis, Gammabis);

