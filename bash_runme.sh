#!/bin/bash

##### find /run/media/smuraru/Tank1/GO/ds_o1/spions/mdps/MG_* -maxdepth 1 -type f -print0 | parallel -0 bash le_script

export GMX_MAXBACKUP=-1
FILES=/run/media/smuraru/Tank1/GO/ds_o1/spions/mdps/MG_*

function doit() {
  ###args
  ###  : @required string f; \
  ##f=/run/media/smuraru/Tank1/GO/ds_o1/spions/mdps/
  chosen=$(ls /run/media/smuraru/Tank1/GO/ds_o1/spions/mdps/ | shuf -n 1)
  move_it=/run/media/smuraru/Tank1/GO/ds_o1/spions/mdps/
  move_it+=$chosen
  mv $move_it /run/media/smuraru/Tank1/GO/ds_o1/spions/mdps/used/
  f=/run/media/smuraru/Tank1/GO/ds_o1/spions/mdps/used/
  f+=$chosen
  echo "Processing $f file..."; \
  n="/run/media/smuraru/Tank1/GO/ds_o1/spions/tprs/"; \
  n+=$(echo $f | sed -r "s/.+\/(.+)\..+/\1/"); \
  n+=".tpr"; \
  xvg_dna_mg="/run/media/smuraru/Tank1/GO/ds_o1/xvg_go_mg/"; \
  xvg_dna_mg+=$(echo $f | sed -r "s/.+\/(.+)\..+/\1/"); \
  xvg_dna_mg+="_MG.xvg"; \
  frame="/run/media/smuraru/Tank1/GO/ds_o1/spions/"; \
  t=$(echo $f | cut -d'_' -f 5); \
  t=$(echo $t | cut -d't' -f 2); \
  t=$(echo $t | cut -d'.' -f 1); \
  t+=".xtc"; \
  frame+=$t; \
  e_out="/run/media/smuraru/Tank1/GO/ds_o1/energies/"; \
  e_out+=$(echo $f | sed -r "s/.+\/(.+)\..+/\1/"); \
  e_out+=".edr"; \
  gmx grompp -f $f -o $n -p topology_local.top -n le_special_ions_MG.ndx -maxwarn 2 -r nvt.gro -c nvt.gro; \
  gmx mdrun -s $n -rerun $frame -e $e_out; \
  echo "63" | gmx energy -f $e_out -o $xvg_dna_mg; \
}

export -f doit
## parallel -j0 doit ::: spions/mdps/*t1*
#parallel -j20 -N 20 --pipe parallel 'doit' ::: {1..60205}

for number in {1..60205}
do
  sem -j +0 "doit"
done
sem --wait

## for f in $FILES
## do
#  sem -j+0 doit $f
##   parallel -j0 doit ::: $f
## done
### sem --wait
## wait
## echo "All done"

#  xvg_go_mg="/run/media/smuraru/Tank1/GO/ds_o1/xvg_go_mg/"
#  xvg_go_mg+=$(echo $f | sed -r "s/.+\/(.+)\..+/\1/")
#  xvg_go_mg+="_MG.xvg"

#  echo $frame

#  echo "63" | gmx energy -f ener.edr -o $xvg_go_mg
