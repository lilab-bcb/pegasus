import unittest

import pegasus as pg
import anndata
import numpy as np
from scipy.sparse import csr_matrix


class TestPreprocessing(unittest.TestCase):
    def test_log_norm(self):
        X = csr_matrix([[1, 11], [2, 20], [5, 6]])
        adata = anndata.AnnData(X)
        adata.var["robust"] = True
        pg.tools.log_norm(adata, 10)
        np.testing.assert_allclose(
            np.expm1(adata.X.toarray()).sum(axis=1), 10, rtol=1e-6, atol=0
        )


if __name__ == "__main__":
    unittest.main()
