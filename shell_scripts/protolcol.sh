#xtc to pdb
system=adna+adna
ref=adna+adna.pdb
input=.xtc
output=adna+adna.pdb
gmx trjconv -s ${ref} -f ${input} -o ${output}

find_pair ${ref} ${system}.inp

x3dna_ensemble analyze -b ${system}.inp -e ${output} -o adna+adna.ensemble.out

x3dna_ensemble extract -f adna+adna.ensemble.out -l

x3dna_ensemble extract -f adna+adna.ensemble.out -p

analyze --abi adna+adna.inp


