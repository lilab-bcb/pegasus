from pegasusio.commands.AggregateMatrix import AggregateMatrix as aggregate_matrix
aggregate_matrix.__doc__ = aggregate_matrix.__doc__.replace('pegasusio', 'pegasus')

from demuxEM.commands.DemuxEM import DemuxEM as demuxEM
demuxEM.__doc__ = demuxEM.__doc__.replace('demuxEM ', 'pegasus demuxEM ')

from .Clustering import Clustering as cluster
from .DeAnalysis import DeAnalysis as de_analysis
from .AnnotateCluster import AnnotateCluster as annotate_cluster
from .Plotting import Plotting as plot
# from .SubClustering import SubClustering as subcluster
from .CheckSampleIndexes import CheckSampleIndexes as check_indexes
from .FindMarkers import FindMarkers as find_markers
from .SCPOutput import SCPOutput as scp_output
