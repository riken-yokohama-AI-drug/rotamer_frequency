# rotamer_frequency
rotamer_frequency is a rotamer energy calculator for 
a set of residue(s) of a protein molecule.

rotamer_frequency provides:
- Rotamer energy calculation of a set of specified residue(s).
- Rotamer energy difference of a set of specified residue(s) between two specified protein structures.
- Classification of rotamer of a set of specified residue(s).

## Requirements
- python 3.9
- pandas 1.5.3
- biopandas 0.2.7
- biopython 1.84
- pymol-open-source 3.0.0
- pyyaml

## How to run rotamer_frequency

 $ python rotamer_frequency.py ChainID:OneLetterAminoAcidResidueNumber,...  pdb_file_path pdb_file_name

### Positional arguments:
  - residue_number, Specify the residue_number
  - model_list, Specify the path to the model list by a text file
  - xray, Specify the xray

### optional arguments:
  - -h, --help            show this help message and exit
  - --chi_atom CHI_ATOM, -c CHI_ATOM, 
                        Specify the chi_atom.
  - --rotamer_lib ROTAMER_LIB, -l ROTAMER_LIB, 
                        Specify the rotamer_lib.
  - --rotameric_def ROTAMERIC_DEF, -r ROTAMERIC_DEF, 
                        Specify the rotamer_dif.
  - --non_rotameric_def NON_ROTAMERIC_DEF, -n NON_ROTAMERIC_DEF, 
                        Specify the rotamer_lib.
  - --is_mutate, -m,  A flag if the input is mutated compated to the reference structure.

### example outputs
  ```text
  $ cat pdb_file_path
  ./4XVS_H.pdb
　$ python ../src/rotamer_frequency.py H:M74,H:F29 pdb_file_path 4XVS_H.pdb
  [RotamerFrequency] Decoy                     : ./4XVS_H.pdb
  [RotamerFrequency] Chain id                  : H
  [RotamerFrequency] Residue no                : 74
  [RotamerFrequency] Residue name(Xray)        : MET
  [RotamerFrequency] Number of chi(Xray)       : 3
  [RotamerFrequency] Number of chi(Mod)        : 3
  [RotamerFrequency] Phi(Xray)                 : -59.3 deg
  [RotamerFrequency] Phi(Mod)                  : -59.3 deg
  [RotamerFrequency] Psi(Xray)                 : -34.9 deg
  [RotamerFrequency] Psi(Mod)                  : -34.9 deg
  [RotamerFrequency] Chi_1 angle_diff          : 0.0
  [RotamerFrequency] Chi_2 angle_diff          : 0.0
  [RotamerFrequency] Chi_3 angle_diff          : 0.0
  [RotamerFrequency] Chi_1 class(Xray)         : 3(g-)(-62.6 deg)
  [RotamerFrequency] Chi_2 class(Xray)         : 3(g-)(-58.7 deg)
  [RotamerFrequency] Chi_3 class(Xray)         : 3(g-)(-65.9 deg)
  [RotamerFrequency] Chi_1 class(Mod)          : 3(g-)(-62.6 deg)
  [RotamerFrequency] Chi_2 class(Mod)          : 3(g-)(-58.7 deg)
  [RotamerFrequency] Chi_3 class(Mod)          : 3(g-)(-65.9 deg)
  [RotamerFrequency] Probability(Xray)         : 0.183205
  [RotamerFrequency] Probability.max(Xray)     : 0.193528
  [RotamerFrequency] Probability(Mod)          : 0.183205
  [RotamerFrequency] Probability.max(Mod)      : 0.193528
  [RotamerFrequency] Chain id                  : H
  [RotamerFrequency] Residue no                : 29
  [RotamerFrequency] Residue name(Xray)        : PHE
  [RotamerFrequency] Number of chi(Xray)       : 2
  [RotamerFrequency] Number of chi(Mod)        : 2
  [RotamerFrequency] Phi(Xray)                 : -48.9 deg
  [RotamerFrequency] Phi(Mod)                  : -48.9 deg
  [RotamerFrequency] Psi(Xray)                 : -39.2 deg
  [RotamerFrequency] Psi(Mod)                  : -39.2 deg
  [RotamerFrequency] Chi_1 angle_diff          : 0.0
  [RotamerFrequency] Chi_2 angle_diff          : 0.0
  [RotamerFrequency] Chi_1 class(Xray)         : 2(t)(-178.3 deg)
  [RotamerFrequency] Chi_2 class(Xray)         : 6(-)(58.2 deg)
  [RotamerFrequency] Chi_1 class(Mod)          : 2(t)(-178.3 deg)
  [RotamerFrequency] Chi_2 class(Mod)          : 6(-)(58.2 deg)
  [RotamerFrequency] Probability(Xray)         : 0.102493
  [RotamerFrequency] Probability.max(Xray)     : 0.457625
  [RotamerFrequency] Probability(Mod)          : 0.102493
  [RotamerFrequency] Probability.max(Mod)      : 0.457625
  [RotamerFrequency] Angle_diff summary        : 0.0
  [RotamerFrequency] Probability total(Xray)   : 0.018777
  [RotamerFrequency] Probability total(Mod)    : 0.018777
  [RotamerFrequency] Rotamer Energy(Xray)      : 0.9 (kcal/mol)
  [RotamerFrequency] Rotamer Energy(Mod)       : 0.9 (kcal/mol)
  [RotamerFrequency] Rotamer Energy_diff(Mod)  : 0.0 (kcal/mol)
  [RotamerFrequency] Probability log(Mod/Xray) : 0.0
  ```
   

## Smooth Backbone-Dependent Rotamer Library 2010 (Shapovalov2011)
- The backbose-dependent rotamer library was
  retrieved from Shapovalov2011 (http://dunbrack.fccc.edu/lab/bbdep2010).
- The rotamer library compiles relationship of rotamer conformations
 and the frequencies in experimental structures.
- The default location of the directory is input/dunbrack2010-library/.

## How rotamer_frequency works
- The rotamer energy is calculated by -R T ln(probability⁄(maximum probability)),
where R and T are the gas constant of 0.001987 kcal/(mol K) and temperature of 300 K
- The maximum probability is selected from probability values defined
 for the identical backbone geometry.
- Rotamer probability of multiple residues is defined by
the joint probability assuming that they are independent.

## Citation

```text
@article{Chiba2024,
  doi = {TBD},
  url = {TBD},
  year = {2024},
  month = TBD,
  publisher = {TBD},
  author = {Shuntaro Chiba and Tsutomu Yamane and Yasushi Okuno and Mitsunori Ikeguchi and Masateru Ohta},
  title = {TBD},
  journal = {TBD}
}
```
## Citation for the Smooth Backbone-Dependent Rotamer Library 2010

```text
@article{Shapovalov2011,
  doi = {https://doi.org/10.1016/j.str.2011.03.019},
  url = {http://dunbrack.fccc.edu/lab/bbdep2010},
  year = {2011},
  month = {March},
  publisher = {Elsevier},
  author = {Maxim V. Shapovalov and Roland L. Dunbrack, Jr.},
  title = {A smoothed backbone-dependent rotamer library for proteins derived from adaptive kernel density estimates and regressions},
  journal = {Structure}
}
```
