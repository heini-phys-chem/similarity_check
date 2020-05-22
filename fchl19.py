#!/usr/bin/env python3

import sys
import random
from datetime import datetime
import numpy as np
from copy import deepcopy
import qml
from qml.kernels import get_global_kernel
from qml.representations import generate_fchl_acsf
from qml.math import cho_solve
from qml.kernels import l2_distance
import itertools
from time import time

# Function to parse datafile to a dictionary
def get_energies(filename):
	""" Returns a dictionary with heats of formation for each xyz-file.
	"""

	f = open(filename, "r")
	lines = f.readlines()
	f.close()

	energies = dict()

	for line in lines:
		tokens = line.split()

		xyz_name = tokens[2]
		hof = float(tokens[3])

#		if hof < 100 and hof > 0: energies[xyz_name] = hof
		energies[xyz_name] = hof

	return energies

if __name__ == "__main__":
	data  = get_energies(sys.argv[1])
	mols = []

	for xyz_file in sorted(data.keys()):
		mol = qml.Compound()
		mol.read_xyz(xyz_file)
		mol.properties = data[xyz_file]
		mol.name = xyz_file
		mols.append(mol)

	x = []
	q = []

	list_of_elements = [1,5,6,7,8,9,17,35]
	for mol in mols:
		x1 = generate_fchl_acsf(mol.nuclear_charges, mol.coordinates, gradients=False, pad=21, elements=list_of_elements)
		x.append(x1)
		q.append(mol.nuclear_charges)

	X		= np.array(x)
	Q		= q

	K = get_global_kernel(X, X, Q, Q, .64)
	Y = np.asarray([ mol.properties for mol in mols ])

	lst = [mol.name for mol in mols]
	lst_old = len(lst)

	for i in range(len(Y)):
		for j in range(i+1,len(Y)):
			if K[i][j] > 0.9999 and (np.abs(Y[i]-Y[j])*627.509) < 1e-4 and mols[j].name in lst:
				lst.remove(mols[j].name)
				print(mols[i].name, mols[j].name, K[i][j], (np.abs(Y[i]-Y[j])*627.509))

	print("{}/{} uniq conformers".format(len(lst), lst_old))


