# hom
Computes the homology of spaces composed of combinations of simpler spaces, using only default Python libraries.

To use it, simply invoke the command `python3 hom.py`.  You will be presented with a prompt asking for homology / cohomology, a description of the space
you'd like to construct, and a description of the coefficients you'd like to use for computing the homology.

Example usage, computing the cohomology of a torus (the product of two circles) with coefficients in `R`:

```
Input H(*) <space description>; coeff:
H* S1 x S1; R
  dim        | 0  1    2 
  cohomology | R  R+R  R
```

Homology and cohomology are both supported; to use either, type `H ` (homology) or `H* ` (cohomology) at the beginning of your input.

To build a space, write an expression building it out of simpler spaces.  Supported base spaces are:
* `Sn` (n-sphere)
* `RPn` (n-dimensional real projective space) 
* `CPn` (n-dimensional complex projective space) 
* `Dn` (n-disks)
* `*` (a point)  

Additionally, you may combine spaces using the operations 
* `x` (product)
* `v` (wedge sum / 1-point union)
* `$` (suspension)
* `_` (disjoint union).

To use coefficients besides `Z`, you may put a semicolon at the end of your prompt and then write a description of the desired coefficients.
The supported Abelian groups / fields are `Z`, `Zn`, `R`, `Q`, `0`, and direct sums (`+`) of these.


The program uses hardcoded versions of Tor, Ext, and Hom functors, as well as a hardcoded tensor product and hardcoded base space homologies.  Thus it may have issues with strange coefficients.

It uses these to apply the Kunneth and UCT theorems from algebraic topology, as well as smaller theorems about how wedge sum / suspension / disjoint union affect homology.

Parentheses and composition are also supported, and the program includes a simple tokenizer and parser to build a tree out of these operations.
