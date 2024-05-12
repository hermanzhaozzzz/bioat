# demo distance_matrix.txt
"""
 , A, B, C, D, E, F
A, 0, 5, 4, 7, 6, 8
B, 5, 0, 7,10, 9,11
C, 4, 7, 0, 7, 6, 8
D, 7,10, 7, 0, 5, 9
E, 6, 9, 6, 5, 0, 8
F, 8,11, 8, 9, 8, 0
"""
# The algorithm currently runs in O( N^3 ) time (empirically, O( N^2.8 ) or so)
# https://github.com/ivo-pontes/NeighborJoiningTree/blob/master/NeighbourJoining.py
