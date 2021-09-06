# Fast Sparse-Matrix Multiplier
Multiply two sparse matrices either in Jagged Diagonal Storage (JDS) or Coordinate Scheme (COOS) format and return the resulting matrix in the same format. COOS is compacter than JDS for storing super-sparse matrices. Otherwise, JDS is the ideal storage format.

The methods used are shown in [Ausarbeitung.pdf](Ausarbeitung.pdf) (German).
Run the following for more information:

`./main --help`

### Usage:
##### COOS
`./main -c file1.txt file2.txt`
##### JDS
`./main -j file1.txt file2.txt`

### Authors:<br>
Abdelrahman A., Omar A., Alexander A.
