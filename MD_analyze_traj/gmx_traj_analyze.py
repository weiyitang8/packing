import sys, os, numpy as np, subprocess as sp



for time in np.linspace(0, 400000,9):

	cmd = 'echo 1 | gmx trjconv -s prod.tpr -f prod_wrap.xtc -o frame_{:03}.pdb -dump {}' .format(int(time/1000), int(time))
	print(cmd)
	sp.call(cmd,shell=1)


