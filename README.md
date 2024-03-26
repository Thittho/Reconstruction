Description
--
This repository contains Magma code for reconstructing non-hyperelliptic genus 4 curves as the locus of a quadratic and cubic form in P^3.
A direct method for the reconstruction of non-hyperelliptic genus 3 curves is also provided.

The file minimization.m contains a routine that reduces the size of non-hyperelliptic curves of genus 4 given by a quadratic and cubic form with integer or rational coefficients. 

Prerequisites
--
You need an installation of Magma.

Installation
--
You can enable the functionality of this package in Magma by attaching the quartic/magma/spec file with AttachSpec. To make this independent of the directory in which you find yourself, and to active this on startup by default, you may want to indicate the relative path in your `~/.magmarc` file, by adding the line
```
AttachSpec("~/Reconstruction/magma/spec");
```

Usage
--
Examples that test the routines in this package can be found in the directory
[`examples`](examples) (a full list of intrinsics is [here](intrinsics.md)).

Basic use of the package is as follows.

* Given invariants of a (generic smooth) non-hyperelliptic genus 4 curve C realized as the intersection of a quadric and a cubic in P^3, return a quadric and a cubic whose intersection is isomorphic to C.

```
QQ := Rationals();
R<x,y,z,w> := PolynomialRing(QQ, 4);
Quadric1 := &+[Random(50)*m : m in MonomialsOfDegree(R,2)];
Cubic1 := &+[Random(50)*m : m in MonomialsOfDegree(R,3)];
J1, Wgt := InvariantsGenus4Curves(Quadric1, Cubic1 : normalize := true);
J1 := WPSMinimize(Wgt, J1);
Quadric2, Cubic2 := ReconstructionGenus4(J1);
IsIsomorphicG4(Quadric2, Cubic2, Quadric1, Cubic1);
```

* Given the Dixmier-Ohno invariants of a (generic smooth) plane quartic C, return a plane quartic which is isomorphic to C.

```
R<x,y,z> := PolynomialRing(QQ, 3);
F := &+[Random(50)*m : m in MonomialsOfDegree(R,4)];
DO, Wgt := DixmierOhnoInvariants(F);
G := ReconstructionGenus3(DO);
IsIsomorphicTernaryQuartics(F,R!G);
```

Citing this work
--
If this code was helpful to your research, please cite the following preprint:

Credits
--
The skeleton of this README.md file was copied from [https://github.com/JRSijsling/quartic/blob/main/README.md)](https://github.com/JRSijsling/quartic/blob/main/README.md).
