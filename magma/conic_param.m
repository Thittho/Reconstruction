/*
 *  Computes the parametrization of a conic, copied from the package https://github.com/JRSijsling/quartic.
 *
 *  Distributed under the terms of the GNU Lesser General Public License (L-GPL)
 *                  http://www.gnu.org/licenses/
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  Copyright:
 *  2020 R. Lercier, C. Ritzenthaler & J.R. Sijsling
 */
 
function ConicParametrization(L : RationalPoint := true, RandomLine := true, Legendre := false, B := 100)
    /* B is the maximal height of the integral coefficients of the intersecting line. */

    K := BaseRing(Parent(L));
    P := ProjectiveSpace(K, 2); x := P.1; y := P.2; z := P.3;
    C := Conic(P, L);

    /* Can we find a rational point on this conic ? */
    if (Type(K) eq FldFin) or (RationalPoint
	and ((Type(K) in {FldRat, FldFin, RngInt}) or
	     (Type(K) eq FldAlg and (
	            Degree(K) eq 1 or IsSimple(K))) or
	     (Type(K) eq FldFunRat and (
		    Type(BaseField(K)) eq FldRat or
		    ISA(Type(BaseField(K)),FldNum) or
		    (IsFinite(BaseField(K)) and Characteristic(BaseField(K)) ne 2 ))) or

	     (ISA(Type(K), FldFunG) and Characteristic(K) ne 2)))
	then
	HasRatPoint, Pt := HasRationalPoint(C);
	if HasRatPoint then
            vprintf Hyperelliptic, 2 : "Conic has a rational point\n";
	    return Parametrization(C, Place(Pt), Curve(ProjectiveSpace(K, 1)));
	end if;
	vprintf Hyperelliptic, 2 : "Conic has no rational point\n";
    end if;

    /* Since we have no rational point on it, let us construct a quadratic extension that contains one */
    if Legendre then
        LC, LMmap := LegendreModel(C); LL := DefiningPolynomial(LC);
    else
        LC := C; LL := DefiningPolynomial(LC);
    end if;

    if RandomLine then
        D := [-B..B];
        repeat
            c1 := Random(D);
            c2 := Random(D);
            c3 := Random(D);
        until c2 ne 0;
        R<t> := PolynomialRing(K);
        h := hom<Parent(LL) -> R | [R.1, -(c1/c2)*R.1 - 1, 1]>;
        S := ExtensionField < K, t | h(LL) >;
        Pt := [ S | S.1, (-c1/c2)*S.1 - 1, 1];
    else
        a := K ! MonomialCoefficient(LL, x^2); b := K ! MonomialCoefficient(LL, x*z); c := K ! MonomialCoefficient(LL, z^2);
        S := ExtensionField < K, x | a*x^2 + b*x + c >;
        Pt := [ S | S.1, 0, 1];
    end if;

    CS := Conic(ProjectiveSpace(S, 2), L);
    if Legendre then
        Pt := [ Evaluate(p, Pt) : p in DefiningPolynomials(Inverse(LMmap)) ];
    end if;
    return Parametrization( CS, Place(CS!Pt), Curve(ProjectiveSpace(S, 1)) );

end function;

