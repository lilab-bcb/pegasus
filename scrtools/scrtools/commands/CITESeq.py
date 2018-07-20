from .Base import Base
from ..pipeline import run_cite_seq_pipeline

class CITESeq(Base):
    """
Run CITE-Seq pipeline.

Usage:
  scrtools cite_seq --max-cells <number> [options] <input_file> <output_name>
  scrtools cite_seq -h

Arguments:
  input_file       Input file in h5at format.
  output_name      Output file name. All outputs will use it as the prefix.

Options:
  -p <number>, --threads <number>                  Number of threads. [default: 1]
  --genome <genome>                                Genome name. [default: GRCh38]
  
  --max-cells <number>                             Only top <number> cells sorted by UMI counts.
  --counts-per-cell-after <number>                 Total counts per cell after normalization. [default: 1e5]
  
  --random-state <seed>                            Random number generator seed. [default: 0]

  --nDC <number>                                   Number of diffusion components. [default: 50]
  --diffmap-alpha <alpha>                          Power parameter for diffusion-based pseudotime. [default: 0.5]
  --diffmap-K <K>                                  Number of neighbors used for constructing affinity matrix. [default: 100]
  --diffmap-full-speed                             For the sake of reproducibility, we only run one thread for building kNN indices. Turn on this option will allow multiple threads to be used for index building. However, it will also reduce reproducibility due to the racing between multiple threads.

  --run-louvain                                    Run louvain clustering algorithm.
  --louvain-resolution <resolution>                Resolution parameter for the louvain clustering algorithm. [default: 1.3]
  --louvain-affinity <affinity>                    Affinity matrix to be used. Could be 'W_norm', 'W_diffmap', or 'W_diffmap_norm'. [default: W_norm]

  --run-approximated-louvain                       Run approximated louvain clustering algorithm.
  --approx-louvain-ninit <number>                  Number of Kmeans tries. [default: 20]
  --approx-louvain-nclusters <number>              Number of clusters for Kmeans initialization. [default: 30]
  --approx-louvain-resolution <resolution>.        Resolution parameter for louvain. [default: 1.3]

  --run-tsne                                       Run multi-core tSNE for visualization.
  --tsne-perplexity <perplexity>                   tSNE's perplexity parameter. [default: 30]

  -h, --help                                       Print out help information.

Outputs:
  output_name.h5ad        Output file in h5ad format. The clustering results are stored in the 'obs' field (e.g. 'louvain_labels' for louvain cluster labels). The PCA, tSNE and diffusion map coordinates are stored in the 'obsm' field.
  
Examples:
  scrtools cite_seq -p 20 --run-louvain --run-tsne manton_bm.h5at manton_bm
    """

    def execute(self):
        kwargs = {
            'n_jobs' : int(self.args['--threads']),
            'genome' : self.args['--genome'],
            'max_cells' : int(self.args['--max-cells']),
            'norm_count' : float(self.args['--counts-per-cell-after']),
            'random_state' : int(self.args['--random-state']),

            'nDC' : int(self.args['--nDC']),
            'diffmap_alpha' : float(self.args['--diffmap-alpha']),
            'diffmap_K' : int(self.args['--diffmap-K']),
            'diffmap_full_speed' : self.args['--diffmap-full-speed'],

            'run_louvain' : self.args['--run-louvain'],
            'louvain_resolution' : float(self.args['--louvain-resolution']),
            'louvain_affinity' : self.args['--louvain-affinity'],

            'run_approx_louvain' : self.args['--run-approximated-louvain'],
            'approx_louvain_ninit' : int(self.args['--approx-louvain-ninit']),
            'approx_louvain_nclusters' : int(self.args['--approx-louvain-nclusters']),
            'approx_louvain_resolution' : float(self.args['--approx-louvain-resolution']),

            'run_tsne' : self.args['--run-tsne'],
            'tsne_perplexity' : float(self.args['--tsne-perplexity'])
        }

        run_cite_seq_pipeline(self.args['<input_file>'], self.args['<output_name>'], **kwargs)
