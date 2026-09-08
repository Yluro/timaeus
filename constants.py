"""Static reference data shared across the SymmetryMeasurements plugin."""

# Element symbols grouped by block. Used to recognize plausible coordination centres.
S_METALS = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra']
P_METALS = ['Al', 'Ga', 'In', 'Tl', 'Ge', 'Sn', 'Pb', 'Sb', 'Bi', 'Po']
D3 = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"]
D4 = ["Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"]
D5 = ["Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"]
LN = ["La", "Ce", "Nd", "Pr", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
AN = ["Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]
METALS = S_METALS + P_METALS + D3 + D4 + D5 + LN + AN

# Reference geometry indices accepted by SHAPE 2.1 for each vertex count,
# taken from the SHAPE 2.1 documentation.
REF_SHAPE_DICT = {2: "1 2 3",
                  3: "1 2 3 4",
                  4: "1 2 3 4",
                  5: "1 2 3 4 5",
                  6: "1 2 3 4 5",
                  7: "1 2 3 4 5 6 7",
                  8: "1 2 3 4 5 6 7 8 9 10 11 12 13",
                  9: "1 2 3 4 5 6 7 8 9 10 11 12 13",
                  10: "1 2 3 4 5 6 7 8 9 10 11 12 13",
                  11: "1 2 3 4 5 6 7",
                  12: "1 2 3 4 5 6 7 8 9 10 11 12 13",
                  20: "1",
                  24: "1 2",
                  48: "1",
                  60: "1"
                  }

shape21_citation = 'Llunell, M., Casanova, D., Cirera, J., Alemany P., Alvarez, S. (2013). SHAPE (Version 2.1). Universitat de Barcelona.'
octadist_citation = 'Ketkaew, R., Tantirungrotechai, Y., Harding, P., Chastanet, G., Guionneau, P., Marchivie, M., & Harding, D. J. (2021). OctaDist: a tool for calculating distortion parameters in spin crossover and coordination complexes. Dalton Transactions, 50(3), 1086-1096.'
