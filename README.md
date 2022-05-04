# O(3P) + HBr -> OH + Br PES

The O(3P) + HBr --> OH + Br MRCI+Q/CBS(aug-cc-pVnZ(-PP); n = Q,5)+SO 3A"
potential surface of A. G. S. de Oliveira-Filho, F. R. Ornellas
and K. A. Peterson

The Reproducing Kernel Hilbert Space method of Ho and Rabitz is used
to interpolate the surface consisting of 1110 geometries spanning
O-H-Br angles of 60-180 deg.  This surface does not contain the
H + BrO arrangement and no claims of accuracy are made for O-H-Br
angles smaller than 60 deg.


## USAGE:

### Input:
      r1 = R(O-H)
      r2 = R(H-Br)
      r3 = R(Br-O).

### Output:
      value = Energy, relative to the asymptotic reactants
              valley (O + HBr(r_e)), in hartree
      dr1 = Derivative of the potential with respect to r1 (R(O-H))
            in hartree/bohr
      dr2 = Derivative of the potential with respect to r2 (R(H-Br))
            in hartree/bohr
      dr3 = Derivative of the potential with respect to r3 (R(O-Br))
            in hartree/bohr

#### NOTE:
    Before any actual potential energy calculations are made, a single
    call to myprepot must be made:
      `call myprepot`

    Later, the potential energy is computed by calling mypot:
      `call mypot(r1,r2,r3,value,dr1,dr2,dr3)`

   The parameters are read from the `o_hbr_a_pp.par` file

   
# References
1. A. G. S. de Oliveira-Filho, F. R. Ornellas, K. A. Peterson, Accurate ab initio potential energy surfaces for the 3A’’ and 3A’ electronic states of the O(3P)+HBr system. J. Chem. Phys. 136, 174316 (2012). (https://doi.org/10.1063/1.4705428)
2. A. G. S. de Oliveira-Filho, F. R. Ornellas, K. A. Peterson, S. L. Mielke, Thermal rate constants for the O(3P) + HBr and O(3P) + DBr reactions: transition-state theory and quantum mechanical calculations. J. Phys. Chem. A 117, 12703–12710 (2013). (https://doi.org/10.1021/jp4090684)

# Geometries, energies, and gradients for testing

## Saddle point
```
      r1                          r2                        r3
 2.6289099999999999        2.8486750000000001        5.1064215904894388
      energy
 7.9866812204093529E-003
      dr1                         dr2                       dr3
 1.5320730578638475E-007   1.5109533511581397E-006  -1.5577335148869720E-007
```

## Reactants side vdW well
```
      r1                          r2                        r3
 4.3746320000000001        2.6855400000000000        7.0601719999996142
      energy
-2.5958471320488563E-003
      dr1                         dr2                       dr3
 1.4650519388234480E-003   1.4648796113941348E-003  -1.4650322069311608E-00
```

## Products  side vdW well
```
      r1                          r2                        r3
 1.8372599999999999        4.8222319999999996        4.4678939513580840
      energy
-3.3549004949264938E-002
      dr1                         dr2                       dr3
 6.5041529057150577E-007  -4.8876198921465885E-010   6.6615884752874166E-009
```
