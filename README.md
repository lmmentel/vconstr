# README file for Vconstr


## basinfo file

The file is written by the subroutine `WRITEBASINFO` located in the `gamess.src` file in the modified
version of GAMESS(US) code for dsfun calculations.

### Format of the file

Name   | Type | Dim    | Fmt           | Description
-------------------------------------------------------------
TITLE  | CHAR |     80 |           a80 | Title
NAT    | INT  |      1 |           i25 | Number of atoms
ICH    | INT  |      1 |           i25 | Charge
MUL    | INT  |      1 |           i25 | Multiplicity
NUM    | INT  |      1 |           i25 | Total number of basis functions 
NX     | INT  |      1 |           i25 | Total number of gaussians 
NE     | INT  |      1 |           i25 | Number of electrons
NA     | INT  |      1 |           i25 | Number of electrons with alpha spin (NE+MUL-1)/2
NB     | INT  |      1 |           i25 | Number of electrons with beta spin (NE-MUL+1)/2
NSHELL | INT  |      1 |           i25 | Number of basis set shells
NPRIMI | INT  |      1 |           i25 | Number of primitive functions in the basis set
ZAN    | REAL |    NAT | '(3(e25.15))' | Nuclear charges for each atom
C      | REAL | 3, NAT | '(3(e25.15))' | Cartesian coordiantes of each atom 
IMIN   | INT  |    NAT |    '(3(i25))' |
IMAX   | INT  |    NAT |    '(3(i25))' |
EVEC   | REAL |      3 | '(3(e25.15))' | External electric fieled components
KATOM  | INT  | NSHELL |    '(3(i25))' | TELLS WHICH ATOM THE SHELL IS CENTERED ON, NORMALLY MORE THAN ONE SHELL EXISTS ON EVERY ATOM.
INTYP  | INT  | NSHELL |    '(3(i25))' | 
EX     | REAL | NPRIMI | '(3(e25.15))' | Gaussian exponents, for every symmetry unique primitive
C1     | REAL | NPRIMI | '(3(e25.15))' | 
C2     | REAL | NPRIMI | '(3(e25.15))' | 
KNG    | INT  | NSHELL |    '(3(i25))' | IS THE NUMBER OF GAUSSIANS IN THIS SHELL.  THEIR DATA ARE STORED CONSECUTIVELY BEGINNING AT THE -KSTART- VALUE.
KLOC   | INT  | NSHELL |    '(3(i25))' | GIVES THE LOCATION OF THIS SHELL IN THE TOTAL AO BASIS, PLEASE READ THE EXAMPLEcin gamess inputa.src file
