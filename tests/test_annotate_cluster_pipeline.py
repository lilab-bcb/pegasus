import unittest
import os

import sccloud.commands


class TestAnnotateClusterPipeline(unittest.TestCase):
    def tearDown(self):
        os.path.exists("annotate_output.txt") and os.remove("annotate_output.txt")

    def test_annotate(self):
        cmd = sccloud.commands.annotate_cluster(
            [
                "annotate_cluster",
                "--marker-file",
                os.path.join("tests", "scCloud-test-data", "input", "markers.json"),
                "--de-test", "fisher",
                os.path.join("tests", "scCloud-test-data", "output", "test_de.h5ad"),
                "annotate_output.txt"
            ]
        )
        cmd.execute()
        with open(
            os.path.join("tests", "scCloud-test-data", "output", "annotate_output.txt")
        ) as f:
            data = f.readlines()
        with open("annotate_output.txt") as f:
            test_data = f.readlines()
        self.assertEqual(data, test_data)


if __name__ == "__main__":
    unittest.main()
