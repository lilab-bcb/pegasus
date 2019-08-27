import unittest
import os
import shutil
import sccloud.commands
from .test_util import assert_excel_equal


class TestDePipeline(unittest.TestCase):
    # def tearDown(self):
    #     os.path.exists("test_de.xlsx") and os.remove("test_de.xlsx")
    #     os.path.exists("test_de.h5ad") and os.remove("test_de.h5ad")

    def test_de_analysis(self):
        # de_analysis modifies h5ad file
        shutil.copy(
            os.path.join("tests", "scCloud-test-data", "output", "test_cluster.h5ad"),
            "test_de.h5ad",
        )
        cmd = sccloud.commands.de_analysis(
            [
                "de_analysis",
                "test_de.h5ad",
                "test_de.xlsx",
                "--fisher",
                "--mwu",
                "--auc",
                "--t",
                "--labels",
                "leiden_labels",
            ]
        )
        cmd.execute()
        assert_excel_equal(
            self,
            os.path.join("tests", "scCloud-test-data", "output", "test_de.xlsx"),
            "test_de.xlsx",
        )


if __name__ == "__main__":
    unittest.main()
