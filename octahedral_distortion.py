"""A reimplementation of the OctaDist distortion parameters, plus tau and mu.

Opposite faces and vertices are identified topologically, from a convex hull that
is expected to be combinatorially equivalent to an octahedron. Structures with
more than three coplanar vertices can degenerate the hull and make that
identification fail.
"""
import os
from itertools import combinations, permutations

import matplotlib.pyplot as plt
import numpy as np
import olx
from scipy.spatial import ConvexHull

from selection import MolecularStructure

# A six-coordinate centre: one central atom and six donor atoms.
EXPECTED_ATOMS = 7


class CalcDistortion:
    def __init__(self, coordination_structure: MolecularStructure):

        if len(coordination_structure) != EXPECTED_ATOMS:
            raise ValueError(f'Expected {EXPECTED_ATOMS} atoms (a centre and six donors), '
                             f'found {len(coordination_structure)}.')

        self.labels = coordination_structure.labels
        self.coords = np.array(coordination_structure.coords, dtype=np.float64)

        self.central_atom = self.coords[0]
        self.vertices = self.coords[1:]

        # Vectors pointing from metal to ligand
        self.vectors = np.array([coord - self.central_atom for coord in self.vertices])

        self.bond_distances = self.calc_bond_distances()
        self.mean_bond_distance = float(np.mean(self.bond_distances))

        # Angles is a sorted array; for an octahedron the last three are the trans angles.
        self.angles = self._all_angles()
        self.cis_angles = self.angles[:-3]
        self.trans_angles = self.angles[-3:]

        self.convex_hull = ConvexHull(self.vertices)
        self.volume = float(self.convex_hull.volume)

        # Array of arrays with indices of self.vertices that make a triangle.
        self.faces = self.convex_hull.simplices
        self.normals = self.calc_normals()

        # Each entry pairs the indices of two faces / two vertices that sit opposite
        # each other: [[[1 2 3] [4 5 6]] ... ]
        self.opposite_faces = self._opposite_faces()
        self.opposite_vertices = self._opposite_vertices()

        # Final calculations
        self.zeta = self.calc_zeta()
        self.delta = self.calc_delta()
        self.sigma = self.calc_sigma()
        self.tau = self.calc_tau()
        self.theta = self.calc_theta()
        self.mu = self.calc_mu()

    def calc_bond_distances(self):
        ds = [np.linalg.norm(v) for v in self.vectors]
        return np.array(ds, dtype=np.float64)

    def _all_angles(self):
        angles = []
        for i, j in combinations(range(len(self.vectors)), 2):
            v1 = self.vectors[i] / np.linalg.norm(self.vectors[i])
            v2 = self.vectors[j] / np.linalg.norm(self.vectors[j])
            angle = np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1, 1)))
            angles.append(angle)

        angles.sort()
        return np.array(angles)

    def calc_zeta(self) -> float:
        deviations = [np.abs(d - self.mean_bond_distance) for d in self.bond_distances]
        return float(np.sum(deviations))

    def calc_delta(self) -> float:
        delta = np.sum(np.power(
            (self.bond_distances - self.mean_bond_distance) / self.mean_bond_distance, 2)) / 6
        return float(delta)

    def calc_sigma(self) -> float:
        return float(np.sum(np.abs(90 - self.cis_angles)))

    def calc_tau(self) -> float:
        return float(np.sum(np.abs(180 - self.trans_angles)))

    def calc_theta(self) -> float:
        thetas = []

        for face_pair in self.opposite_faces:
            angles = []
            # Assign front and back faces
            front, back = face_pair
            # The normal is calculated from the front face
            normal = self.normals[front]

            # For each index of the front face
            for i in self.faces[front]:
                # Get the vector i and its projection onto the front face
                v_i = self.vectors[i]
                p_i = self._project_onto_plane(v_i, normal)
                p_i = p_i / np.linalg.norm(p_i)

                for j in self.faces[back]:
                    pair = [i, j]
                    pair_r = [j, i]
                    is_opposite = any(np.array_equal(pair, ov) or np.array_equal(pair_r, ov)
                                      for ov in self.opposite_vertices)
                    if is_opposite:
                        continue

                    v_j = self.vectors[j]
                    p_j = self._project_onto_plane(v_j, normal)
                    p_j = p_j / np.linalg.norm(p_j)
                    angle = np.degrees(np.arccos(np.clip(np.dot(p_i, p_j), -1, 1)))
                    angles.append(np.abs(angle))

            thetas.append(np.sum([np.abs(60 - angle) for angle in angles]))

        return float(np.average(thetas) * 4)

    def calc_mu(self) -> float:
        centroid = self.coords[1:].mean(axis=0)
        return float(np.linalg.norm(centroid - self.coords[0]))

    def calc_normals(self):
        normals = []
        for face in self.faces:
            # A face is a tuple of indices that form a triangle in self.vertices
            p1, p2, p3 = self.vertices[face[0]], self.vertices[face[1]], self.vertices[face[2]]
            # Calculate two in-plane vectors of the plane
            v1 = p2 - p1
            v2 = p3 - p1
            normal = np.cross(v1, v2)
            normal = normal / np.linalg.norm(normal)
            normals.append(normal)
        return np.array(normals)

    def _opposite_faces(self):
        pairs = []
        # For each ordered pair of faces:
        for i, j in permutations(range(len(self.faces)), 2):
            # An empty intersection means the faces share no vertex,
            # so in an octahedron they are opposite each other.
            if not set(self.faces[i]) & set(self.faces[j]):
                pairs.append((i, j))
        return np.array(pairs)

    def _opposite_vertices(self):
        pairs = []
        # For every combination of vertex i and j
        for i, j in combinations(range(len(self.vertices)), 2):
            # Get all faces that contain vertex i
            faces_with_i = [face for face in self.faces if i in face]

            # For each of those faces, get all vertices that form them
            vertices_near_i = set(v for face in faces_with_i for v in face)

            # If j is not in the set, then i and j share no common face,
            # so i and j are opposite.
            if j not in vertices_near_i:
                pairs.append((i, j))

        return np.array(pairs)

    def _project_onto_plane(self, vector, normal):
        return vector - np.dot(vector, normal) * normal

    def print_results(self, file_name=''):
        print('\n' + '=' * 70)
        print(f'Octahedral distortion parameters calculated for {self.labels[0]} in {file_name}')
        print('-' * 70)
        print(f"{'Mean d(M-X)':<12}{self.mean_bond_distance:>12.4f}   {'Ang':<12}")
        print(f"{'Zeta':<12}{self.zeta:>12.4f}   {'Ang':<12}")
        print(f"{'Delta':<12}{self.delta:>12.6f}   {'':<12}")
        print(f"{'Sigma':<12}{self.sigma:>12.2f}   {'deg':<12}")
        print(f"{'Theta':<12}{self.theta:>12.2f}   {'deg':<12}")
        print(f"{'Volume':<12}{self.volume:>12.4f}   {'Ang^3':<12}")
        print('-' * 70)
        print(f"{'Tau':<12}{self.tau:>12.2f}   {'deg':<12}")
        print(f"{'Mu':<12}{self.mu:>12.4f}   {'Ang':<12}")
        print('=' * 70)

    def draw_octahedron(self):
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.tight_layout()

        ax.scatter(self.central_atom[0], self.central_atom[1], self.central_atom[2],
                   color='blue', s=50, zorder=5)
        ax.scatter(self.vertices[:, 0], self.vertices[:, 1], self.vertices[:, 2],
                   color='red', s=50, zorder=5)

        colors = plt.cm.tab10(np.linspace(0, 1, len(self.faces)))

        for i, face in enumerate(self.faces):
            triangle = [self.vertices[face[0]], self.vertices[face[1]], self.vertices[face[2]]]
            poly = Poly3DCollection([triangle], alpha=0.3, facecolor=colors[i], edgecolor='black')
            ax.add_collection3d(poly)

        for label, coord in zip(self.labels, self.coords):
            ax.text(coord[0], coord[1], coord[2], label, ha='right', va='top', zorder=10)

        # Olex2 uses the non-GUI version of matplotlib, so the plot can only be
        # saved and viewed later - there is no interactive 3D view.
        try:
            plt.show()
        except UserWarning:
            plt.close()

        try:
            save_dir = os.path.join(olx.FilePath(), 'Oh_distortion')
            os.makedirs(save_dir, exist_ok=True)
            plt.savefig(os.path.join(save_dir, 'octahedron.png'))
            print(f'Octahedron graph saved to {save_dir}.')
        except Exception as e:
            print(f'Could not save the octahedron graph: {e}')
        finally:
            plt.close()
