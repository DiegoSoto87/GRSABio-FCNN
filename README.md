# GRSABio
Golden Ratio Simulated Annealing Bio-inspired (GRSABio) Copyright (C) 2025 Dr. Hernán Peraza Vazquez and Dr. Diego Arturo Soto Monterrubio,
Golden Ratio Simulated annealing (GRSA) Copyright (C) 2020  Dr. Juan Paulo Sánchez Hernández, and Dr. Juan Frausto Solis
Copyright (C) 2005 Frank Eisenmenger, U.H.E. Hansmann, Shura Hayryan, Chin-Ku Hu This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.

README
# GRSABio-FCNN
The Golden Ratio Simulated Annealing Bio-inspired (GRSABio) algorithm is a hybrid approach that combines the Golden Ratio Simulated Annealing method with the Jumping Spider Optimization Algorithm (JSOA) to refine the structures of initial peptide models assembled from fragments predicted by a Convolutional Neural Network (FCNN). The fragment-assembled templates used for testing are provided in the EXAMPLES folder, which contains 60 peptide templates for evaluating the algorithm's performance.
_________________________________________________________________________________________________________
1. Installation.
To install the GRSABio algorithm you need fortran 90 and an operating system with linux. GRSABio algorithms come per folder for each one which contains the SMMP package, inside the each folder there is a makefile needed you compile the SMMP.

In the linux terminal you enter the file path, for example cd home/user/GRSABio-FCNN. Then in the terminal type make clean, this instruction cleans the object files. Next in terminal type make, this create the new object files.

For each GRSABio-FCNN algorithm you need to do the make clean and make, in order to install the program.

2. Configuration.
The EXAMPLES folder contains the sequence file (1a13.seq) and template file (1a13.var) for each algorithm, this files are necessary to run the program with the instance 1a13. The instance parameters are in grsabio.f

Run algorithm.
For run the algorithm program you need type ./smmp in the linux terminal. example: ~/GRSABio-FCNN$ ./smmp

The program will start to run and show the energy with its respective temperature until the execution is finished.

Output files
Threedimensional structure: performs the prediction of the three-dimensional structure in a *.pdb format

The Energy of the structure: contains the energy of prediction of the three-dimensional structure in *.txt format

Troubleshooting For any problem you can contact diego_060787@hotmail.com
The GRSABio package is This software is the Golden Ratio Simulated Annealing Bio-inspired (GRSABio). It is based on the SMMP package. The GRSABio software is based on the Golden Ratio simulated annealing and Jumping Spider Optimization Algorithm algorithms described in the next papers:

Peraza-Vázquez, H.; Peña-Delgado, A.; Ranjan, P.; Barde, C.; Choubey, A.; Morales-Cepeda, A. B. A bio-inspired method for845
mathematical optimization inspired by Arachnida salticidade. Mathematics, 2021, 10, 102.

Juan Frausto-Solis,Juan Paulo Sánchez-Hernández, Mishael Sánchez-Pérez, Ernesto Liñan García, Golden Ratio Simulated Annealing for Protein Folding Problem. International Journal of Compututational Methods , Vol 12, No. 6,1550037, 2015. Doi: https://doi.org/10.1142/S0219876215500371.

The following lines are from the SMMP README file:

The SMMP package is designed for molecular simulation of linear proteins within the standard geometry model. The package contains various modern Monte Carlo algorithms and energy minimization routines. The energy of the protein can be calculated by exploiting one of three force fields: ECEPP/2, ECEPP/3 or FLEX. Two subroutines are included for approximating protein-solvent interaction by means of calculating the solvent accessible area of atomic groups (either analytically or in an approximate way). More information on SMMP can be found in the manual (file `manual.ps'), our papers (F. Eisenmenger, U.H.E. Hansmann, S. Hayryan and C.K. Hu, [SMMP] A Modern Package for Protein Simulations, Comp.Phys. Comm. (2001), 138 (2001) 192-212; [SMMP 1.1] - An Enhanced Open-Source Software Package, submitted), and
on the SMMP homepage: http://www.smmp05.net

SMMP is offered as open code. We encourage users to re-write code, or add their own modules, whenever they see the need. For terms of use of SMMP, please see the manual (manual.pdf) and the license agreement (license.txt). By viewing or using any part of SMMP you agree with the terms of the license. Any suggestions for improvement of the code or reports on bugs are welcome.

Please send your remarks to:

-Frank Eisenmenger: eisenmenger@fmp-berlin.de
-Ulrich H.E. Hansmann: hansmann@mtu.edu
-Shura Hayryan: shura@phys.sinica.edu.tw
-Chin-Kun Hu: huck@phys.sinica.edu.tw
-Jan H. Meinke j.meinke@fz-juelich.de

INSTALLATION: After uncompressing and un-tar-ing the SMMP package into a separate directory, the following files should be in that directory:

*README This README-file

*SMMP A sub-directory containing:

-lib.sh2 Library file with ECEPP/2 parameters
-lib.sh3 Library file with ECEPP/3 parameters
-lib.flex Library file with FLEX parameters
-charges File with charges (needed only for FLEX)
-tes.dat File with tesselation points needed for approximating the solvent accessible surface area of atoms

*doc A sub-directory containing:

-manual.pdf The manual in PDF-format describing the details of various subroutines and installation of SMMP
-license.txt Ascii-file with the license agreement
-angle_defs.svg Scalable vector graphic of definition of angles
-atom_numbering.svg Scalable vector graphic of numbering of atoms
-dihedral_defs.svg Scalable vector graphic of definition of dihedrals
-Makefile Makefile to build documentation
-manual.lyx Original source of manual
-manual.tex Lyx file converted to tex suitable for pdflatex

*Makefile Produces the executable file

*INCL.H Defines global parameters and common-blocks

*INCP.H Defines parameters and common-blocks necessary for reading PDB-structures

*incl_lund.h Defines parameters and common blocks used with the Lund force field and BGS

*init_molecule.f Wrapping function that constructs the start configuration of a molecule

*init_energy.f Wrapping function that initializes energy parameter

*redseq.f Reads sequence of amino acids (e.g. enkefa.seq)

*redvar.f Reads initial configuration of molecule (e.g. enkefa.var)

*getmol.f Assembles data from libraries

*bldmol.f Builds up the protein atomic coordinates

*addend.f Adds end-groups

*setmvs.f Determines sets of moving atoms for given variables

*mklist.f Compiles lists of interaction partners

*setvar.f Resets variables and rebuilds the molecule

*dihedr.f Returns dihedral angles and valence angles

*difang.f Calculates the difference and the sum of two angles

*nursvr.f Calculates the residue index of a dihedral angle or atom

*redstr.f String input routines

*pdbread.f Reads amino acid sequence and rotein atomic coordinates from a PDB-file, calculates dihedral angles from PDB coordinates, and builds index field that relates PDB-atoms to SMMP atoms

*energy.f Wrapping function that returns the energy of a current protein configuration

*enyflx.f Calculates internal energy of molecule with FLEX dataset

*enyshe.f Calculates internal energy of molecule with ECEPP datasets

*enyshe_p.f Parallel version of enyshe.f

*enylun.f Calculates energy of molecule using the Lund force field

*enysol.f Calculates solvation energy of molecule using solvent accessible area method (fast, but approximate calculation)

*enysol_p.f Parallel version of enysol.f

*esolan.f Calculates solvation energy of molecule using solvent accessible area method (analytic, but slow calculation)

*eninteract.f Calculates interaction term for multi-molecule simulations based on ECEPP.

*enyreg.f Calculates a constraint energy needed for regularizing PDB-structures

*eyabgs.f Calculates correction term introduced by R. A. Abagyan et al.

*gradient.f Wrapping function that returns the energy gradient vs. dihedral angles for a protein configuration

*opeflx.f Calculates internal energy and partial derivatives vs. dihedral angles for FLEX dataset

*opeshe.f Calculates internal energy and partial derivatives vs. dihedral angles for ECEPP datasets

*opesol.f Calculates analytically the partial derivatives vs. dihedral angles of the solvation energy

*opereg.f Calculates the partial derivatives vs. diheadral angles of the constraint energy term during regularization

*main.f Main program

*regul.f Regularization of PDB-structure into SMMP geometry

*anneal.f For simulated annealing run

*canon.f For canonical Monte Carlo run

*minim.f For minimization of protein potential energy

*mulcan_par.f Calculates multicanonical weights

*mulcan_sim.f For multicanonical simulation run

*mulcan_par_mod.f90 Combines the previous two files in a Fortran module.

*partem_s.f For parallel tempering run on a single-processor machine

*partem_p.f For parallel tempering run on a multiple-processor machine

*main_p.f Replaces `main.f' on multiple-processor machine in parallel tempering runs.

*main_bgl_p.f Version of main_p.f optimized for use with IBM BlueGene/L

*metropolis.f Performs Metropolis updates

*bgs.f Biased Gaussian step.

*minqsn.f Minimization by quasi-Newton method using BFGS-formula

*mincjg.f Minimization by conjugate gradient method

*outvar.f Output of the current conformation (dihedral angles)

*contacts.f Calculates van der Waals contacts

*cnteny.f Calculates atomic contact energy and prints bad contacts

*hbond.f Calculates number of hydrogen bonds in a configuration

*helix.f Measures the number of residues which are part of helix or beta-sheet

*outpdb.f Output of current configuartion in PDB-format

*rgyr.f Measures the radius-of-gyration and end-to-end distance in molecule

*zimmer.f Expresses given configuration in Zimmerman code

*rmsdfun.f Calculates root-mean-square deviation between current SMMP configuration and a reference structure.

*twister.f Mersenne-Twister random number generator

*utilties.f Somer helper function for simulations on multiple processors.

*smmp.pyf Interface of the Python bindings

*universe.py Python package that provides access to global properties of the system.

*protein.py Python package representing atoms, amino acids, and proteins

*algorithms.py Some basic algorithms using the Python bindings.

*main.py Main program in Python

*rmexclpoint.py Utility script to build the Python binding

*restoreexclpoint.py Utility script to build the Python binding

*temperatures Sample temperature file with 32 temperatures.

*energy.f Wrapping function that returns the energy of a current protein configuration

*enyflx.f Calculates internal energy of molecule with FLEX dataset

*enyshe.f Calculates internal energy of molecule with ECEPP datasets

*enyshe_p.f Parallel version of enyshe.f

*enylun.f Calculates energy of molecule using the Lund force field

*enysol.f Calculates solvation energy of molecule using solvent accessible area method (fast, but approximate calculation)

*enysol_p.f Parallel version of enysol.f

*esolan.f Calculates solvation energy of molecule using solvent accessible area method (analytic, but slow calculation)

*eninteract.f Calculates interaction term for multi-molecule simulations based on ECEPP.

*enyreg.f Calculates a constraint energy needed for regularizing PDB-structures

*eyabgs.f Calculates correction term introduced by R. A. Abagyan et al.

*gradient.f Wrapping function that returns the energy gradient vs. dihedral angles for a protein configuration

*opeflx.f Calculates internal energy and partial derivatives vs. dihedral angles for FLEX dataset

*opeshe.f Calculates internal energy and partial derivatives vs. dihedral angles for ECEPP datasets

*opesol.f Calculates analytically the partial derivatives vs. dihedral angles of the solvation energy

*opereg.f Calculates the partial derivatives vs. diheadral angles of the constraint energy term during regularization

*main.f Main program

*regul.f Regularization of PDB-structure into SMMP geometry

*anneal.f For simulated annealing run

*canon.f For canonical Monte Carlo run

*minim.f For minimization of protein potential energy

*mulcan_par.f Calculates multicanonical weights

*mulcan_sim.f For multicanonical simulation run

*mulcan_par_mod.f90 Combines the previous two files in a Fortran module.

*partem_s.f For parallel tempering run on a single-processor machine

*partem_p.f For parallel tempering run on a multiple-processor machine

*main_p.f Replaces `main.f' on multiple-processor machine in parallel tempering runs.

*main_bgl_p.f Version of main_p.f optimized for use with IBM BlueGene/L

*metropolis.f Performs Metropolis updates

*bgs.f Biased Gaussian step.

*minqsn.f Minimization by quasi-Newton method using BFGS-formula

*mincjg.f Minimization by conjugate gradient method

*outvar.f Output of the current conformation (dihedral angles)

*contacts.f Calculates van der Waals contacts

*cnteny.f Calculates atomic contact energy and prints bad contacts

*hbond.f Calculates number of hydrogen bonds in a configuration

*helix.f Measures the number of residues which are part of helix or beta-sheet

*outpdb.f Output of current configuartion in PDB-format

*rgyr.f Measures the radius-of-gyration and end-to-end distance in molecule

*zimmer.f Expresses given configuration in Zimmerman code

*rmsdfun.f Calculates root-mean-square deviation between current SMMP configuration and a reference structure.

*twister.f Mersenne-Twister random number generator

*utilties.f Somer helper function for simulations on multiple processors.

*smmp.pyf Interface of the Python bindings

*universe.py Python package that provides access to global properties of the system.

*protein.py Python package representing atoms, amino acids, and proteins

*algorithms.py Some basic algorithms using the Python bindings.

*main.py Main program in Python

*rmexclpoint.py Utility script to build the Python binding

*restoreexclpoint.py Utility script to build the Python binding

*temperatures Sample temperature file with 32 temperatures.

*EXAMPLES A sub-directory containing:

-enkefa.seq Example sequence file

-enkefa.var Example configuration file

-enkefa.ann Example configuration file for simulated annealing run

-enkefa.ref Example contact matrix file

-abeta.seq Sequence file for Abeta_16-22

-abeta.var Configuration file with global coordinates

-abeta.ref Dummy contact matrix

-1bdd.pdb PDB file of protein A for regularization.

-1bdd.ref Contact matrix of protein A

-1bdd.seq Sequence file of protein A

-1bdd.var Sample configuration of protein A

-1a13.var Sample configuration of protein

-1a13.seq Sequence file of protein

-temperatures Example temperature file for parallel_tempering_s

-temperatures_abeta Example temperature file for parallel_tempering

-smmp.cmd Example shell-script to run smmp

-Makefile Makefile for building the examples.

-annealing.f Example for running simulated annealing

-minimization.f Example for a minimization

-multicanonical.f Example for calculating parameters for a multi canonical simulation.

-parallel_tempering_p.f Example for parallel tempering on a cluster

-parallel_tempering_s.f Example for parallel tempering on a single node

-partem_p.f Parallel tempering routine used with parallel_tempering_p.f The output is different from the default implementation.

-regularization.f Example for regularizing protein A (1bdd)

-Python A sub-directory containing some examples that use the Python bindings.

-annealing.py Python version of the annealing example

-minimization.py Python version of the minimization example

-muca.py Python version for calculating the parameters for a multi-canonical simulation.

-regularization.py Python version of regularization example

-gui_example.py Example for building a graphical user interface for SMMP

-best.pml A PyMol script for rendering best.pdb

-scripts A sub-directory containing some useful scripts:

  +README Short description of the scripts
  +atomprops.py Lists the properties of all atoms in a protein
  +var2pdb.py Takes a sequence and a var file and builds a PDB from it.
  
The whole SMMP package is written in standard FORTRAN language. We have been exploiting it under pgif, ifort, gfortran, and xlf. It should be possible to compile the code with any contemporary fortran compilers. There are no machine dependent routines included in SMMP.

Common blocks and limiting parameters are gathered in special files INCL.H''  and INCP.H'' which are attached to the modules through an include' statement.  In order to install SMMP the user needs to edit the Makefile' and specify the compiler and compiler options which he will use. Executing the make' command will finish installation of SMMP. For compiling the parallel version use make parallel'. The command make doc' build the documentation and make pybind' build the Python bindings if f2py is installed.

  *HOW TO RUN SMMP: SMMP does not include an interpretor of user defined commands. The preparation of a simulation must be done in the 'main' module After     changes of SMMP has to be re-compiled.
  
Alternatively, you can use the Python bindings, which allow for interactive simulations.

The residues which can be used with each parameter set are described in files `lib.sh2', 'lib.sh3' and 'lib.flx', respectively. The file 'charges' is needed for N- and C-terminal residues with FLEX parameters. The directory with these 3 files should be given in string 'libdir', which is assigned in module 'main' or 'pmain'.

SMMP requires as input a file that specifies the sequence of residues in a protein. This sequence can be read either from a PDB-file (the standard format in which protein structures are deposited in the Protein Data Bank) or from a special sequence file. If the sequence is read from a sequence file its first line must start with a '#' and may (or may not) contain the name for the molecule. The residues in the following lines should be named as in the libraries (not case- sensitive). Residue names should be separated by at least one space. An example file ("enkefa.seq") is provided in the sub-directory EXAMPLES.

The initial values for internal variables, i.e. dihedral angles for single bonds, can be calculated either from the atomic coordinates of the PDB-file, or (if the sequence is read from a sequence file) may be provided in a SECOND INPUT file. If this file is not given (or the name of a non-existing file is entered), all variables retain their values given in the libraries. The example file ("enkefa.var") which is provided in the subdirectory EXAMPLES demonstrates the syntax that has to be used:

               residue(s) : variable(s) : value

In the first field the RESIDUE is selected through an {\it INTEGER} number which marks the position of that RESIDUE in the amino acid sequence. The second field lists a string with the name of the VARIABLE, i.e. names the specific dihedral angle. The last field lists the value (a {\it REAL} number) for the VARIABLE and is mandatory. Missing fields are interpreted as 'for all'. Spaces are not significant and are ignored. Empty lines or or lines containing '#' are ignored. A line containing '&' assigns FIXED variable(s), i.e. they will be set to the given value, but will NOT be varied during subsequent changes of the protein configuration.

The following steps summarize how to run SMMP:

*Assign to the character variable 'libdir' the path to the directory containing the standard amino acid residue libraries and the file 'charges'.

*Select the force field and solvation model by setting the four 'sh2', 'epsd' and 'itysol', 'ientyp' to their appropriate values:

  ientyp : 0 => ECEPP2 or ECEPP3 depending on the value of sh2 1 => FLEX 2 => Lund force field 3 => ECEPP with Abagyan corrections
  sh2 =.TRUE. : ECEPP/2 potential, sh2=.FALSE.: ECEPP/3 (Note that the variable `flex' has to be set to .FALSE.)
  epsd=.TRUE. : Distant dependent epsilon(r) epsd=.FALSE.: epsilon=2
  itysol = 0 : Gas phase;
  itysol > 0 : approximation of protein-solvent interactions by means of a solvent accessible surface area approach with stochastic estimation of the accessible   area. itysol < 0 : same as above, but the accessible area is calculated analytically (much SLOWER than itysol > 0).
*Choose a N-terminal and C-terminal group by setting 'grpn' and 'grpc' to approbriate values.

*Choose how the initial input is read in:

  iabin = 0 : read from PDB-file
  iabin = 1 : read from sequence (and configuration) file

*Enter the names of the corresponding file(s). In the example version of the 'main' module this is done through interactive dialog but the user can easily just assign the corresponding names to the character variables 'seqfil' and 'varfil' in the subroutine 'init_molecule'.

*At this point the program is ready for calling the simulation subroutines. In the provided version the energy minimization subroutine is called through 'call minim'. A detailed description of this and other simulation subroutines can be found in the manual (file manual.ps). Normally the simulation subroutines write data in output files, but one can also put output routines such as `outpdb' in 'main'. The minimal output (written into standard output) is the name of the sequence file (extension .seq), name of configuration file (extension .var), and for each residue a list of dihedral angles together with their initial values.

*For parallel tempering jobs on on a multiprocessor system one has to replace 'main' by 'main_p'. The above protocol still applies.

LIMITATIONS: All parameters which limit the usage of SMMP are stored in the file ``INCL.H''. The most important ones are listed below.

mxml=10 max. number of molecules mxrs=500 max. total number of residues mxat=10000 max. total number of atoms mxbd=3 max. number of bonds to following atoms mxvr=mxrs5 max. number of local variables mxms=mxvr3 max. total number of moving sets mxvw=mxat4 max. number of vwd domains mx14=mxat4 max. number of '1-4' partners mxath = 100 max. number of atoms in help-arrays mxvrh=mxath max. number of variables in help mxtyat=18 max. number of energetic atom-types mxhbdo=4 max. types of Hydrogens as donors in HB mxhbac=6 max. types of atoms as acceptors in HB mxtyto=19 max. number of types of torsional potentials nrsty=35 max. number of residue types mxtysol=9 the number of solvation parameters sets

Note also the following restrictions in the current version of SMMP:

  A single amino acid residue can not be simulated with FLEX potential.
  A protein must not start with a prolin residue.

EXAMPLE: Proper installation of SMMP can be tested by running the following example. After compilation of the program package (with make'  command using the provided Makefile') and runing SMMP by typing ./smmp, SMMP will minimize the ECEPP/3 energy of the Met-enkephalin configuration in `EXAMPLES/enkefa.var'.
_____________________________________________________________________________
|NOTE: If the program doesn't start but only shows the error message | |"Killed", you probably don't have enough memory available. You can | |reduce the memory requirement by setting lower limits for mxrs | |(line 10) and mxat (line 11) in the file INCL.H. Setting mxrs=10 and | |mxat=1300 will still run all the examples. |
Running the program leads to the following output:

------------------------------------------------------------ enkefa.out

file with SEQUENCE: ./EXAMPLES/enkefa.seq

file with VARIABLES: ./EXAMPLES/enkefa.var

redvar> Met-Enkephalin: residue 1 Tyr : x1 set -172.590 redvar> Met-Enkephalin: residue 1 Tyr : x2 set 78.710 redvar> Met-Enkephalin: residue 1 Tyr : x6 set -165.880 redvar> Met-Enkephalin: residue 1 Tyr : phi set -86.240 redvar> Met-Enkephalin: residue 2 Gly : psi set 156.180 redvar> Met-Enkephalin: residue 2 Gly : omg set -180.000 redvar> Met-Enkephalin: residue 2 Gly : phi set -154.530 redvar> Met-Enkephalin: residue 3 Gly : psi set 83.640 redvar> Met-Enkephalin: residue 3 Gly : omg set 180.000 redvar> Met-Enkephalin: residue 3 Gly : phi set 83.660 redvar> Met-Enkephalin: residue 4 Phe : psi set -73.860 redvar> Met-Enkephalin: residue 4 Phe : omg set -180.000 redvar> Met-Enkephalin: residue 4 Phe : x1 set 58.790 redvar> Met-Enkephalin: residue 4 Phe : x2 set 94.600 redvar> Met-Enkephalin: residue 4 Phe : phi set -137.040 redvar> Met-Enkephalin: residue 5 Met : psi set 19.330 redvar> Met-Enkephalin: residue 5 Met : omg set -180.000 redvar> Met-Enkephalin: residue 5 Met : x1 set 52.760 redvar> Met-Enkephalin: residue 5 Met : x2 set 175.280 redvar> Met-Enkephalin: residue 5 Met : x3 set -179.830 redvar> Met-Enkephalin: residue 5 Met : x4 set -58.570 redvar> Met-Enkephalin: residue 5 Met : phi set -163.630 redvar> Met-Enkephalin: residue 5 Met : pst set 160.450 redvar> Met-Enkephalin: residue 5 Met : omt set 180.000

Energy BEFORE minimization:

Total: 0.10354E+04 Coulomb: 0.1991E+02 Lennard-Jones: 0.1142E+03 HB: 0.9008E+03 Variables: 0.4511E+00 Solvatation: 0.0000E+00

Step 1: energy 0.489976E+05 ( 0.867243E+13 ) Step 2: energy 0.108241E+03 ( 0.964364E+07 ) Step 3: energy -0.632755E+00 ( 0.498256E+04 ) Step 4: energy -0.739539E+01 ( 0.109979E+05 ) Step 5: energy 0.517957E+03 ( 0.678892E+09 ) Step 6: energy -0.767774E+01 ( 0.143880E+05 ) Step 7: energy -0.818570E+01 ( 0.100364E+05 ) Step 8: energy -0.938401E+01 ( 0.172508E+04 ) Step 9: energy -0.109799E+02 ( 0.105718E+04 ) Step 10: energy -0.116244E+02 ( 0.580821E+03 ) Step 11: energy -0.116969E+02 ( 0.145676E+04 ) Step 12: energy -0.118569E+02 ( 0.856882E+03 ) Step 13: energy -0.119576E+02 ( 0.521381E+03 ) Step 14: energy -0.120897E+02 ( 0.398602E+03 ) Step 15: energy -0.121476E+02 ( 0.222201E+03 ) Step 16: energy -0.121804E+02 ( 0.242613E+03 ) Step 17: energy -0.121975E+02 ( 0.237765E+03 ) Step 18: energy -0.122092E+02 ( 0.217435E+03 ) Step 19: energy -0.122177E+02 ( 0.198576E+03 ) Step 20: energy -0.122259E+02 ( 0.181726E+03 ) Step 21: energy -0.122298E+02 ( 0.178796E+03 ) Step 22: energy -0.122318E+02 ( 0.178968E+03 ) Step 23: energy -0.122352E+02 ( 0.181503E+03 ) Step 24: energy -0.122407E+02 ( 0.184272E+03 ) Step 25: energy -0.122499E+02 ( 0.186336E+03 ) Step 26: energy -0.122634E+02 ( 0.183760E+03 ) Step 27: energy -0.122776E+02 ( 0.171223E+03 ) Step 28: energy -0.122929E+02 ( 0.148718E+03 ) Step 29: energy -0.123122E+02 ( 0.119336E+03 ) Step 30: energy -0.123348E+02 ( 0.847677E+02 ) Step 31: energy -0.123572E+02 ( 0.434069E+02 ) Step 32: energy -0.123733E+02 ( 0.176053E+02 ) Step 33: energy -0.123858E+02 ( 0.941361E+01 ) Step 34: energy -0.123956E+02 ( 0.739923E+01 ) Step 35: energy -0.124014E+02 ( 0.513793E+01 ) Step 36: energy -0.124049E+02 ( 0.385630E+01 ) Step 37: energy -0.124083E+02 ( 0.416421E+01 ) Step 38: energy -0.124122E+02 ( 0.530410E+01 ) Step 39: energy -0.124159E+02 ( 0.493920E+01 ) Step 40: energy -0.124184E+02 ( 0.306992E+01 ) Step 41: energy -0.124201E+02 ( 0.174747E+01 ) Step 42: energy -0.124216E+02 ( 0.110719E+01 ) Step 43: energy -0.124231E+02 ( 0.793785E+00 ) Step 44: energy -0.124242E+02 ( 0.585080E+00 ) Step 45: energy -0.124249E+02 ( 0.454481E+00 ) Step 46: energy -0.124257E+02 ( 0.491193E+00 ) Step 47: energy -0.124268E+02 ( 0.498392E+00 ) Step 48: energy -0.124277E+02 ( 0.297774E+00 ) Step 49: energy -0.124281E+02 ( 0.987143E-01 ) Step 50: energy -0.124282E+02 ( 0.562089E-01 ) Step 51: energy -0.124283E+02 ( 0.472681E-01 ) Step 52: energy -0.124283E+02 ( 0.360471E-01 ) Step 53: energy -0.124284E+02 ( 0.149691E-01 ) Step 54: energy -0.124285E+02 ( 0.376934E-02 ) Step 55: energy -0.124285E+02 ( 0.243872E-02 ) Step 56: energy -0.124285E+02 ( 0.253572E-02 ) Step 57: energy -0.124285E+02 ( 0.176076E-02 ) Step 58: energy -0.124285E+02 ( 0.952400E-03 ) Step 59: energy -0.124285E+02 ( 0.889573E-03 ) Step 60: energy -0.124285E+02 ( 0.747233E-03 ) Step 61: energy -0.124285E+02 ( 0.417036E-03 ) Step 62: energy -0.124285E+02 ( 0.307630E-03 ) Step 63: energy -0.124285E+02 ( 0.474529E-03 ) Step 64: energy -0.124285E+02 ( 0.506211E-03 ) Step 65: energy -0.124285E+02 ( 0.298789E-03 ) Step 66: energy -0.124285E+02 ( 0.125304E-03 ) Step 67: energy -0.124285E+02 ( 0.153436E-03 ) Step 68: energy -0.124285E+02 ( 0.132386E-03 ) Step 69: energy -0.124285E+02 ( 0.661655E-04 ) Step 70: energy -0.124285E+02 ( 0.194869E-04 ) Step 71: energy -0.124285E+02 ( 0.120524E-04 ) Step 72: energy -0.124285E+02 ( 0.641640E-05 ) Step 73: energy -0.124285E+02 ( 0.182141E-05 ) Step 74: energy -0.124285E+02 ( 0.219308E-05 ) Step 75: energy -0.124285E+02 ( 0.348349E-05 ) Step 76: energy -0.124285E+02 ( 0.166901E-05 ) Step 77: energy -0.124285E+02 ( 0.177102E-06 ) Step 78: energy -0.124285E+02 ( 0.663522E-08 ) Step 79: energy -0.124285E+02 ( 0.210693E-08 ) Step 80: energy -0.124285E+02 ( 0.370905E-08 ) Step 81: energy -0.124285E+02 ( 0.324300E-08 ) Step 82: energy -0.124285E+02 ( 0.777363E-09 ) Step 83: energy -0.124285E+02 ( 0.443359E-10 ) Step 84: energy -0.124285E+02 ( 0.322695E-11 ) Step 85: energy -0.124285E+02 ( 0.443359E-10 ) ---- CONVERGENCE ----

Final energies __________________________________________________

Total: -0.12429E+02 Coulomb: 0.2143E+02 Lennard-Jones: -0.2923E+02 HB: -0.6706E+01 Variables: 0.2084E+01 Solvatation: 0.0000E+00

Variables _________________

x1 1 -173.2 ( 0.6) x2 1 79.3 ( 0.6) x6 1 -166.3 ( 0.5) phi 1 -83.1 ( 3.2) psi 2 155.8 ( 0.4) omg 2 -177.1 ( 2.9) phi 2 -154.2 ( 0.3) psi 3 85.8 ( 2.2) omg 3 168.5 ( 11.5) phi 3 83.0 ( 0.7) psi 4 -75.0 ( 1.2) omg 4 -170.0 ( 10.0) x1 4 58.9 ( 0.1) x2 4 94.5 ( 0.1) phi 4 -136.8 ( 0.2) psi 5 19.1 ( 0.2) omg 5 -174.1 ( 5.9) x1 5 52.9 ( 0.1) x2 5 175.3 ( 0.0) x3 5 -179.9 ( 0.0) x4 5 -58.6 ( 0.0) phi 5 -163.4 ( 0.2) pst 5 160.8 ( 0.3) omt 5 -179.8 ( 0.2)

Gradient ______________________________________________________________ 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000


