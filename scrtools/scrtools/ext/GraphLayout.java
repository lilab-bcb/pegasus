import org.gephi.graph.api.Graph;
import org.gephi.graph.api.GraphController;
import org.gephi.graph.api.GraphModel;
import org.gephi.graph.api.Node;
import org.gephi.io.exporter.api.ExportController;
import org.gephi.io.importer.api.Container;
import org.gephi.io.importer.api.EdgeDirectionDefault;
import org.gephi.io.importer.api.ImportController;
import org.gephi.io.processor.plugin.DefaultProcessor;
import org.gephi.layout.plugin.forceAtlas2.ForceAtlas2;
import org.gephi.project.api.ProjectController;
import org.gephi.project.api.Workspace;
import org.openide.util.Lookup;
import org.gephi.layout.plugin.fruchterman.FruchtermanReingold;
import org.gephi.layout.plugin.forceAtlas2.ForceAtlas2;
import org.gephi.layout.plugin.openord.OpenOrdLayout;
import org.gephi.layout.spi.Layout;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;


public class GraphLayout {

	public static void main(String[] args) throws IOException {
		File file = new File(args[0]);
		String output = args[1];
		String layoutString = args[2];
		int nsteps = Integer.parseInt(args[3]);
		int njobs = Integer.parseInt(args[4]); // Add an argument for number of jobs
		ProjectController pc = Lookup.getDefault().lookup(ProjectController.class);

		pc.newProject();
		Workspace workspace = pc.getCurrentWorkspace();
		ImportController importController = Lookup.getDefault().lookup(ImportController.class);
		GraphModel graphModel = Lookup.getDefault().lookup(GraphController.class).getGraphModel();

		Container container = importController.importFile(file);
		container.getLoader().setEdgeDefault(EdgeDirectionDefault.UNDIRECTED);
		Graph g = graphModel.getUndirectedGraph();
		importController.process(container, new DefaultProcessor(), workspace);

		if (layoutString.equals("fa")) {
			ForceAtlas2 fa = new ForceAtlas2(null);
			fa.setGraphModel(graphModel);
			fa.resetPropertiesValues();
			fa.setBarnesHutOptimize(true);
			fa.setThreadsCount(njobs);
			fa.initAlgo();
			for (int i = 0; i < nsteps && fa.canAlgo(); i++) {
				fa.goAlgo();
				if ((i + 1) % 500 == 0) System.out.println("Finished " + String.valueOf(i + 1) + " iterations.");
			}
			fa.endAlgo();
		} else if (layoutString.equals("oo")) {
			OpenOrdLayout oo = new OpenOrdLayout(null);
			oo.setGraphModel(graphModel);
			oo.resetPropertiesValues();
			oo.setNumThreads(njobs);
			oo.setNumIterations(nsteps);
			oo.initAlgo();
			oo.goAlgo();
			oo.endAlgo();
		} else {
			System.err.println("Unknown layout");
		}

		ExportController ec = Lookup.getDefault().lookup(ExportController.class);
		PrintWriter pw = new PrintWriter(new FileWriter(output));
		pw.print("id\tx\ty\n");
		for (Node n : g.getNodes()) {
			pw.print(n.getId());
			pw.print("\t");
			pw.print(n.x());
			pw.print("\t");
			pw.print(n.y());
			pw.print("\n");
		}
		pw.close();
	}
}
