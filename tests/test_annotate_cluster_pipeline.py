import unittest
import os

import sccloud.commands


class TestAnnotateClusterPipeline(unittest.TestCase):
    def tearDown(self):
        os.path.exists("test_annotate.txt") and os.remove("test_annotate.txt")

    def test_annotate(self):
        cmd = sccloud.commands.annotate_cluster(
            [
                "annotate_cluster",
                os.path.join("tests", "scCloud-test-data", "output", "test_de.h5ad"),
                "--annotation",
                os.path.join("tests", "scCloud-test-data", "input", "markers.json"),
            ]
        )
        cmd.execute()
        with open(
            os.path.join("tests", "scCloud-test-data", "output", "test_annotate.txt")
        ) as f:
            data = f.readlines()
        with open("test_annotate.txt") as f:
            test_data = f.readlines()
        self.assertEqual(data, test_data)


if __name__ == "__main__":
    unittest.main()
