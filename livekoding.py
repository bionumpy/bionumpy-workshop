import numpy as np

sequences = np.array([
    [0, 1, 2, 3],
    [10, 11, 12, 13],
    [1, 1, 1, 1]
])


mask = sequences < 10

print(np.sum(mask, axis=0))



