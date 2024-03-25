/////////////////////////////////////////////////////////////////////
////////////* Random example with quadric of rank 4 *////////////////
/////////////////////////////////////////////////////////////////////

R<X,Y,Z,T> := PolynomialRing(Rationals(), 4);

Q1 := &+[Random(1,10)*m : m in [X^2,Y^2,Z^2,T^2]];
E1 := &+[Random(-5,5)*m : m in MonomialsOfDegree(R, 3)];

inv := InvariantsGenus4Curves(Q1, E1);

Q2, E2 := ReconstructionGenus4(inv);

IsIsomorphicGenus4(Q1,E1,Q2,E2);


/////////////////////////////////////////////////////////////////////
////////////* Random example with quadric of rank 3 *////////////////
/////////////////////////////////////////////////////////////////////

R<X,Y,Z,T> := PolynomialRing(Rationals(), 4);

Q1 := &+[Random(-5,5)*m : m in [X^2,X*Y,X*Z,Y^2,Y*Z,Z^2]];
E1 := &+[Random(1,10)*m : m in MonomialsOfDegree(R, 3)];

inv := InvariantsGenus4Curves(Q1, E1);

Q2, E2 := ReconstructionGenus4(inv);
R2 := Parent(Q2);

IsIsomorphicGenus4(Q1,E1,R2!Q2,R2!E2);

