"""Helper classes and functions for dealing with points and distances
in the Riemannian symmetric space SL(n, R) / SO(n).

Eventually, some of the functionality in this module might get folded
into the main geometry_tools package, but for now it is provided
separately.

"""

import numpy as np

from geometry_tools import utils

class SymmetricPoint:
    """Class to represent a point in the symmetric space SO(n,R) / SO(n).
    """
    def __init__(self, matrix=None, a=None, k=None):
        """Specify a point in the symmetric space as either an n*n positive
        definite symmetric matrix, or an orthogonal matrix together
        with a vector of positive eigenvalues.

        If matrix is specified, then k,a cannot be, and vice versa.

        Parameters
        ----------
        matrix : ndarray
            Array of shape (..., n, n) giving an ndarray of symmetric
            matrices, each specifying a point in the symmetric space.

        k : ndarray
            Array of shape (..., n, n) giving an ndarray of orthogonal
            matrices. Each orthogonal matrix k is part of the
            definition of a symmetric matrix of the form k * D * k^T,
            where D is diagonal.

        a : ndarray
            Array of shape (..., n) giving an ndarray of eigenvalue
            vectors for diagonal matrices. Each diagonal matrix is
            part of the definition of a symmetric matrix of the form k
            * D * k^T
        """
        if matrix is None and (k is None or a is None):
            raise ValueError("A matrix or its diagonalization must be specified to initialize a point")

        if (k is None) ^ (a is None):
            raise ValueError("Both K and A must be specified to initialize a point via diagonalization")

        self.matrix = matrix
        self.k = k
        self.a = a

        if a is None or k is None:
            self._compute_coset(matrix)

        if matrix is None:
            self._compute_matrix(a, k)

    def _compute_matrix(self, a, k):
        # compute symmetric matrix from diagonal and orthogonal matrix
        self.matrix = (k @ utils.construct_diagonal(a) @
                       np.swapaxes(k, -1, -2))

    def _compute_coset(self, matrices):
        # diagonalize a symmetric matrix to get an orthogonal and diagonal matrix.

        # WARNING: we take an abs here, which is a sign that round-off
        # error might be affecting the code badly.
        a, self.k = np.linalg.eigh(matrices)
        self.a = np.abs(a)

    def __getitem__(self, val):
        return SymmetricPoint(self.matrix[val], self.a[val], self.k[val])

    def to_origin(self):
        """Compute an isometry taking this point to the "origin"
        (i.e. identity coset) in the symmetric space.

        """
        return (self.k @
         np.sqrt(np.abs(utils.construct_diagonal(1 / self.a))) @
         np.swapaxes(self.k, -1, -2))

    def origin_midpoint(self):
        """Compute a point halfway along the unique geodesic between this
           point and the "origin"
        """
        return SymmetricPoint(a=np.sqrt(self.a), k=self.k)

    def distance_to_origin(self):
        """Compute the distance from this point to the "origin"
        """
        return np.linalg.norm(np.log(np.abs(self.a)), axis=-1)

    def invert_by_transvection(self):
        """Compute an "opposite" point on the other side of the "origin."

        More precisely: if c:(-infty, infty) -> X is the unique
        unit-speed geodesic in X such that c(0) is the "origin" and
        c(t) is self, then compute c(-t).

        """
        d = self.a.shape[-1]
        reverse = utils.permutation_matrix(list(range(d - 1, -1, -1)))
        return SymmetricPoint(k=self.k @ reverse, a=np.flip(1 / self.a, axis=-1))

    @staticmethod
    def from_translation(matrix):
        """Compute a point by translating the "origin" by some element of
        SL(n, R).
        """
        return SymmetricPoint(matrix @ np.swapaxes(matrix, -1, -2))

def apply_isometry(isometry, point, broadcast="pairwise"):
    """Apply an isometry (or several isometries) to a point (or several
    points) in the symmetric space.

    Parameters
    ----------
    isometry : ndarray
        array of shape (..., n, n) giving an ndarray of elements of
        SL(n, R) to apply

    point : SymmetricPoint
        point (or points) to apply the isometry (or isometries)
        to. This should always be a single instance of the
        SymmetricPoint class; the underlying data of that instance
        could be either a single point or many points.

    broadcast : {"pairwise", "elementwise", "pairwise_reversed"}

        Rule for applying isometry to point.

        If "pairwise" (the default), apply each isometry to each
        point, and return a point whose underlying data combines the
        shapes of isometry with (the underlying data of) point.

        If "pairwise_reversed", do the same, but the underlying data
        of the resulting point is combined in the opposite order.

        If "elementwise", apply isometries to points
        element-by-element. If isometry and point have the same shape,
        then the resulting point object has data with this shape also;
        otherwise, it must be possible to broadcast the shapes of
        isometry and point to each other.

    """
    iso_dims = len(isometry.shape) - 2
    pt_dims = len(point.matrix.shape) - 2

    if broadcast == "elementwise":
        result = isometry @ point.matrix @ np.swapaxes(isometry, -1, -2)

    if broadcast == "pairwise_reversed":
        left = utils.matrix_product(isometry, point.matrix, broadcast="pairwise_reversed")
        result = left @ np.swapaxes(isometry, -1, -2)

    if broadcast == "pairwise":
        left = utils.matrix_product(isometry, point.matrix, broadcast="pairwise")
        right = np.expand_dims(np.swapaxes(isometry, -1, -2),
                               axis=tuple(range(iso_dims, iso_dims + pt_dims)))
        result = left @ right

    return SymmetricPoint(result)

def riemannian_distance(p1, p2, broadcast="pairwise"):
    """Compute Riemannian distance between a pair of points.

    Parameters
    ----------
    p1, p2 : SymmetricPoint
        Points to compute distances between

    broadcast : {"pairwise", "elementwise", "pairwise_reversed"} Rule
        for computing distances between arrays of points; see the
        docstring for apply_isometry.

    """
    g = p1.to_origin()

    g_translate = apply_isometry(g, p2, broadcast=broadcast)

    singular_values, _ = np.linalg.eigh(g_translate)
    vector_valued_dist = np.log(np.abs(singular_values))
    return np.linalg.norm(vector_valued_dist, axis=-1)

def zeta_angle_from_origin(p1, p2, zeta, broadcast="pairwise"):
    """Compute 'zeta-angles' at the origin distance between a pair of points.

    Whenever p is a regular point, and zeta is a regular tangent
    vector at the origin, there is a unique tangent vector at o of the
    form z_p = k_p * zeta, that is contained in a Weyl sector based at
    the origin and containing p. When p1, p2 are regular points in the
    symmetric space, this function computes the angle between the
    tangent vectors z_p1 and z_p2.

    Parameters
    ----------
    p1, p2 : SymmetricPoint
        Points to compute zeta-angles for.

    zeta : ndarray
        diagonal matrix with entries summing to zero, specifying a
        direction in the (standard) Cartan subalgebra of SL(n, R)

    broadcast : {"pairwise", "elementwise", "pairwise_reversed"} Rule
        for computing zeta-angles between arrays of points; see the
        docstring for apply_isometry.

    """
    zeta_trace = np.trace(zeta @ zeta)

    p1_zeta = p1.k @ zeta @ np.swapaxes(p1.k, -1, -2)
    p2_zeta = p2.k @ zeta @ np.swapaxes(p2.k, -1, -2)

    product = utils.matrix_product(p1_zeta, p2_zeta,
                                   broadcast=broadcast)

    return np.trace(product, axis1=-1, axis2=-2) / zeta_trace
