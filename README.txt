program runs with commands
--> python main.py sequences.txt global blosum62.txt -10 -5
--> python main.py sequences.txt global pam.txt -10 -5
or
--> python main.py sequences.txt local blosum62.txt -10 -5
--> python main.py sequences.txt local pam.txt -10 -5

in general

--> python file.py sequence_file algorithm_type scoringmatrix_file gap gap_extension

sequences.txt --> file of sequences
local/global --> algorithm type
blosum62/pam --> scoring matrix 
gap --> gap penalty
gap extension --> gap extension penalty 

scoring matrix format should be same as blosum62.txt 
if format is different, program can be crashed