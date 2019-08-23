import unittest
import sccloud.commands
import os
from .test_util import assert_excel_equal


class TestFindMarkersPipeline(unittest.TestCase):
    def tearDown(self):
        os.path.exists("test.markers.xlsx") and os.remove("test.markers.xlsx")

    def test_find_markers(self):
        command = sccloud.commands.find_markers(
            [
                "find_markers",
                os.path.join("tests", "scCloud-test-data", "output", "test_de.h5ad"),
                "test.markers.xlsx",
                "--labels",
                "leiden_labels",
                "--remove-ribo",
            ]
        )
        command.execute()
        assert_excel_equal(
            self,
            os.path.join("tests", "scCloud-test-data", "output", "test.markers.xlsx"),
            "test.markers.xlsx",
        )


if __name__ == "__main__":
    unittest.main()
