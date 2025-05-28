from .Base import Base
from pegasus.pipeline import run_qc_pipeline


class QC(Base):
    """
Calculate QC metrics.

Usage:
    pegasus qc [options] <input_file> <output_name>
    pegasus qc -h

Arguments:
    input_file       Input file in either 'zarr', 'h5ad', 'loom', '10x', 'mtx', 'csv', 'tsv' or 'fcs' format.
    output_name      Output file name. All outputs will use it as the prefix.

Options:
    -p <number>, --threads <number>                  Number of threads. [default: 1]
    --genome <genome_name>                           If sample count matrix is in either DGE, mtx, csv, tsv or loom format, use <genome_name> as the genome reference name.
    --min-genes <number>                             Only keep cells with at least <number> of genes.
    --max-genes <number>                             Only keep cells with less than <number> of genes.
    --min-counts <number>                            Only keep cells with at least <number> of UMIs.
    --max-counts <number>                            Only keep cells with less than <number> of UMIs.
    --percent-mito <percent>                         Only keep cells with mitochondrial percent less than <percent>%.
    --gene-percent-cells <percent>                   Only use genes that are expressed in at least <percent>% of cells to select variable genes. [default: 0.05]

    -h, --help                                       Print out help information.

Outputs:
    output_name.h5ad                     Output file in Zarr format.

Examples:
    pegasus qc -p 4 manton_bm_10x.h5 manton_bm
    """

    def execute(self):
        kwargs = {
            "genome": self.args["--genome"],
            "min_genes": int(self.args["--min-genes"]) if self.args["--min-genes"] is not None else None,
            "max_genes": int(self.args["--max-genes"]) if self.args["--max-genes"] is not None else None,
            "min_umis": int(self.args["--min-counts"]) if self.args["--min-counts"] is not None else None,
            "max_umis": int(self.args["--max-counts"]) if self.args["--max-counts"] is not None else None,
            "percent_mito": float(self.args["--percent-mito"]) if self.args["--percent-mito"] is not None else None,
            "gene_percent_cells": float(self.args["--gene-percent-cells"]),
        }

        import logging
        logger = logging.getLogger("pegasus")
        fh = logging.FileHandler(f"{self.args['<output_name>']}.log")
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        fh.setFormatter(formatter)
        logger.addHandler(fh)

        run_qc_pipeline(self.args["<input_file>"], self.args["<output_name>"], **kwargs)
