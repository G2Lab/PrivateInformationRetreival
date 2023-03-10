# PrivateInformationRetrival
Private Information Retrieval and Analysis of Genomic Data with Homomorphic Encryption

## Building Project

Download and compile the latest stable version of HELib (make sure you install the [patchelf](https://github.com/NixOS/patchelf) and [m4](https://www.gnu.org/software/m4/) dependency as well).
Use the guide [here](https://github.com/homenc/HElib/blob/master/INSTALL.md) and follow **Option 1**

If you followed the instructions and created a `mylibs` folder in the root directory, run `./scripts/make.sh` and then `make` to compile our repository.

## Sample Run

We have include a sample DB and query script to demonstrate the functionalities of this project. After running the make command, run `./bin/main` to see our sample output (which will be the same as below).

```
Initialising context object...
m = 17293, p = 131, phi(m) = 17292
  ord(p) = 3
  normBnd = 1.27324
  polyNormBnd = 1.27324
  factors = [17293]
  generator 2 has order (== Z_m^*) of 5764
r = 1
nslots = 5764
hwt = 0
ctxtPrimes = [6,7,8,9,10,11,12,13]
specialPrimes = [14,15,16]
number of bits = 595

security level = 92.0971

Security: 92.0971
Num slots: 5764
Printing DB
-----------------------------------------------------
|snp1|snp2|ALS|
--------------
|0    0   0   |
|0    1   0   |
|0    2   0   |
|1    0   0   |
|1    1   1   |
|1    2   1   |
|2    0   0   |
|2    1   1   |
|2    2   0   |
|0    1   0   |
Running sample queries:
-----------------------------------------------------
Running Counting query (snp 0 = 0 and snp 1 = 1)
Count: 2
Running Counting query (snp 0 = 0 and snp 1 = 2)
Count: 1
Running Counting query (snp 0 = 0 or snp 1 = 1)
Count: 6
Running MAF query filter (snp 0 = 0 or snp 1 = 1), target snp = 0
Nom: 3
Dom: 12
Computed MAF: 0.25
Running Distribution query (alpha: snp 0 = 2 and snp 1 = 5)
0, 5, 10, 2, 7, 12, 4, 9, 14, 5
```

