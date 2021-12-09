#!/bin/bash
#update 10/12/2021


python rewrite_ca_restr_OVERWRITE.py
python rewrite_posres.py


gmx grompp -f minim.mdp -c system_solv_ions.gro  -r system_solv_ions.gro -p topol.top  -o em.tpr
gmx mdrun -v -deffnm em 
gmx make_ndx -f em.gro < save.txt



gmx grompp -f equilibrationNVT.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -n index.ndx   
gmx mdrun -deffnm nvt -s nvt.tpr -v      -nt 8 -ntmpi 4 -nb gpu -bonded cpu -pme gpu -npme 1 -nstlist 200

gmx grompp -f equilibrationNPT.mdp -c nvt.gro -r em.gro -p topol.top -o npt.tpr -n index.ndx
gmx mdrun -deffnm npt -s npt.tpr -v     -nt 8 -ntmpi 4 -nb gpu -bonded cpu -pme gpu -npme 1 -nstlist 200

gmx grompp -f production.mdp -c npt.gro -p topol.top -o prod.tpr -n index.ndx
gmx mdrun -deffnm prod -s prod.tpr  -v    -nt 8 -ntmpi 4 -nb gpu -bonded cpu -pme gpu -npme 1 -nstlist 200 #200ns first then 400ns or 500ns


gmx trjconv -s prod.tpr -f prod.xtc -o prod_wrap.xtc -pbc mol -ur compact -center < wrap_text.txt
gmx trjconv -s prod.tpr -f prod_wrap.xtc -o lastFrame_testrun.pdb -dump 200000 < wrap_text.txt



gmx convert-tpr -s prod.tpr -extend 200000 -o prod.tpr
gmx mdrun -v -deffnm prod -s prod.tpr -cpi prod.cpt -nt 8 -ntmpi 4 -nb gpu -bonded cpu -pme gpu -npme 1 -nstlist 200 #extend 200ns



gmx trjconv -s prod.tpr -f prod.xtc -o prod_wrap.xtc -pbc mol -ur compact -center < wrap_text.txt
gmx trjconv -s prod.tpr -f prod_wrap.xtc -o lastFrame_run.pdb -dump 400000 < wrap_text.txt
