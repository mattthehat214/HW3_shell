import math


class MatrixSolver: #just grouping my methods in a class for organization
    def __init__(self, A, b):
        self.A = A
        self.b = b

    def is_symmetric(self):
        """Check if the matrix A is symmetric."""
        for i in range(len(self.A)):
            for j in range(i + 1, len(self.A)):
                if self.A[i][j] != self.A[j][i]:
                    return False
        return True

    def is_positive_definite(self):
        """Check if the matrix A is positive definite."""
        for i in range(1, len(self.A) + 1):
            if self.determinant([row[:i] for row in self.A[:i]]) <= 0:
                return False
        return True

    def determinant(self, matrix):
        """Calculate the determinant of a matrix."""
        if len(matrix) == 1:
            return matrix[0][0]
        if len(matrix) == 2:
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

        det = 0
        for c in range(len(matrix)):
            det += ((-1) ** c) * matrix[0][c] * self.determinant(self.minor(matrix, 0, c))
        return det

    def minor(self, matrix, i, j):
        """Calculate the minor of a matrix element."""
        return [row[:j] + row[j + 1:] for row in (matrix[:i] + matrix[i + 1:])]

    def cholesky_decompose(self):
        """Perform a Cholesky decomposition on A."""
        n = len(self.A)
        L = [[0.0] * n for _ in range(n)]

        for i in range(n):
            for j in range(i + 1):
                sum_k = sum(L[i][k] * L[j][k] for k in range(j))

                if i == j:  # Diagonal elements
                    L[i][j] = math.sqrt(max(self.A[i][i] - sum_k, 0))
                else:
                    L[i][j] = (self.A[i][j] - sum_k) / L[j][j]
        return L

    def lu_decompose(self):
        """Perform LU decomposition using Doolittle's method."""
        n = len(self.A)
        L = [[0.0] * n for _ in range(n)]
        U = [[0.0] * n for _ in range(n)]

        for i in range(n):
            L[i][i] = 1.0

            for j in range(i, n):  # Calculate U
                sum_upper = sum(L[i][k] * U[k][j] for k in range(i))
                U[i][j] = self.A[i][j] - sum_upper

            for j in range(i + 1, n):  # Calculate L
                sum_lower = sum(L[j][k] * U[k][i] for k in range(i))
                L[j][i] = (self.A[j][i] - sum_lower) / U[i][i]

        return L, U

    def forward_substitution(self, L, b):
        """Perform forward substitution to solve Ly = b."""
        y = [0.0] * len(b)
        for i in range(len(b)):
            y[i] = (b[i] - sum(L[i][j] * y[j] for j in range(i))) / L[i][i]
        return y

    def backward_substitution(self, U, y):
        """Perform backward substitution to solve Ux = y."""
        x = [0.0] * len(y)
        for i in reversed(range(len(y))):
            x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, len(y)))) / U[i][i]
        return x

    def solve(self):
        """Solve the matrix equation Ax = b using either Cholesky or LU decomposition."""
        if self.is_symmetric() and self.is_positive_definite():
            L = self.cholesky_decompose()
            y = self.forward_substitution(L, self.b)
            x = self.backward_substitution(L, y)  # This should likely be backward_substitution(L, y) as before
            method = 'Cholesky'
        else:
            L, U = self.lu_decompose()
            y = self.forward_substitution(L, self.b)
            x = self.backward_substitution(U, y)
            method = 'LU'

        return x, method


# Define the matrix A and vector b for the first system of equations
A1 = [
    [1, -1, 3, 2],
    [-1, 5, -5, -2],
    [3, -5, 19, 3],
    [2, -2, 3, 21]
]
b1 = [15, -35, 94, 1]

# Define the matrix A and vector b for the second system of equations
A2 = [
    [4, 2, 4, 0],
    [2, 2, 3, 2],
    [4, 3, 6, 3],
    [0, 2, 3, 9]
]
b2 = [20, 36, 60, 122]


# Create instances of the MatrixSolver for each system
solver1 = MatrixSolver(A1, b1)
solver2 = MatrixSolver(A2, b2)

# Solve the first system
solution1, method_used1 = solver1.solve()
print(f"Solution for the first system using {method_used1} method: {solution1}")

# Solve the second system
solution2, method_used2 = solver2.solve()
print(f"Solution for the second system using {method_used2} method: {solution2}")
