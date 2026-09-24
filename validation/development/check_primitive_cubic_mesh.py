"""Count cubic operations preserving a shifted primitive-BCC MP mesh."""

import itertools

import numpy as np


def main():
    primitive = np.array([[1, 0, 0.5], [0, 1, 0.5], [0, 0, 0.5]])
    points = set(itertools.product([1, 3], repeat=3))
    preserved = 0
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product([-1, 1], repeat=3):
            cartesian = np.eye(3, dtype=int)[:, permutation] @ np.diag(signs)
            fractional = np.linalg.solve(primitive, cartesian @ primitive)
            rotation = np.rint(fractional).astype(int)
            assert np.max(np.abs(fractional - rotation)) == 0
            preserved += all(
                tuple((rotation.T @ k) % 4) in points for k in points
            )
    assert preserved == 8
    print(f"Primitive BCC shifted 2x2x2 MP: {preserved}/48 operations preserve the set")
    print("Gamma: 48/48 operations preserve the set")


if __name__ == "__main__":
    main()
