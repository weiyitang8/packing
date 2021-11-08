

import collections
import math

import numpy as np

from Bio.PDB.kdtrees import KDTree

__all__ = ["ShrakeRupley_nsc"]

_ENTITY_HIERARCHY = {
    "A": 0,
    "R": 1,
    "C": 2,
    "M": 3,
    "S": 4,
}

ATOMIC_RADII = collections.defaultdict(lambda x: 2.0)
ATOMIC_RADII.update(
    {
        "H": 1.200,
        "HE": 1.400,
        "C": 1.700,
        "N": 1.550,
        "O": 1.520,
        "F": 1.470,
        "NA": 2.270,
        "MG": 1.730,
        "P": 1.800,
        "S": 1.800,
        "CL": 1.750,
        "K": 2.750,
        "CA": 2.310,
        "NI": 1.630,
        "CU": 1.400,
        "ZN": 1.390,
        "SE": 1.900,
        "BR": 1.850,
        "CD": 1.580,
        "I": 1.980,
        "HG": 1.550,
    }
)


class ShrakeRupley_nsc:


    def __init__(self, probe_radius=1.40, n_points=100, radii_dict=None):

        self.probe_radius = float(probe_radius)
        self.n_points = n_points
        # Update radii list with user provided lists.
        self.radii_dict = ATOMIC_RADII.copy()
        if radii_dict is not None:
            self.radii_dict.update(radii_dict)

        # Pre-compute reference sphere
        self._sphere = self._compute_sphere()

    def _compute_sphere(self):
        """Return the 3D coordinates of n points on a sphere.
        Uses the golden spiral algorithm to place points 'evenly' on the sphere
        surface. We compute this once and then move the sphere to the centroid
        of each atom as we compute the ASAs.
        """
        n = self.n_points

        dl = np.pi * (3 - 5 ** 0.5)
        dz = 2.0 / n

        longitude = 0
        z = 1 - dz / 2

        coords = np.zeros((n, 3), dtype=np.float32)
        for k in range(n):
            r = (1 - z * z) ** 0.5
            coords[k, 0] = math.cos(longitude) * r
            coords[k, 1] = math.sin(longitude) * r
            coords[k, 2] = z
            z -= dz
            longitude += dl

        return coords #coords = np.zeros((n, 3), dtype=np.float32)

    def compute(self, atoms):


        n_atoms = len(atoms)


        # Get coordinates as a numpy array
        # We trust DisorderedAtom and friends to pick representatives.
        coords = np.array([a.coord for a in atoms], dtype=np.float64)

        # Pre-compute atom neighbors using KDTree
        kdt = KDTree(coords, 10)

        # Pre-compute radius * probe table
        radii_dict = self.radii_dict
        radii = np.array([radii_dict[a.element] for a in atoms], dtype=np.float64)
        radii += self.probe_radius
        twice_maxradii = np.max(radii) * 2

        # Calculate ASAs
        asa_array = np.zeros((n_atoms, 1), dtype=np.int)
        ptset = set(range(self.n_points))
        for i in range(n_atoms):

            r_i = radii[i]

            # Move sphere to atom
            s_on_i = (np.array(self._sphere, copy=True) * r_i) + coords[i]
            available_set = ptset.copy()

            # KDtree for sphere points
            kdt_sphere = KDTree(s_on_i, 10)

            # Iterate over neighbors of atom i
            for jj in kdt.search(coords[i], twice_maxradii):
                j = jj.index
                if i == j:
                    continue

                if jj.radius < (r_i + radii[j]):
                    # Remove overlapping points on sphere from available set
                    available_set -= {
                        pt.index for pt in kdt_sphere.search(coords[j], radii[j])
                    }

            asa_array[i] = len(available_set)  # update counts

        # Convert accessible point count to surface area in A**2
        f = radii * radii * (4 * np.pi / self.n_points)
        asa_array = asa_array * f[:, np.newaxis]

        # Set atom .sasa
        for i, atom in enumerate(atoms):
            atom.sasa = asa_array[i, 0]


















                