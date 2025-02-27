from scipy.spatial import cKDTree

pos = [[0, 0], [100, 100], [10, 10]]
kdtree = cKDTree(pos)
pos = [[0, 0], [1, 1], [10, 10]]
print(kdtree.query(pos[0], k=2))
