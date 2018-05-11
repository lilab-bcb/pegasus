workflow scrtools_graph_layout {
    # in h5ad format
	File input_file
	String output_name
	# fa (ForceAtlas2), fr (FruchtermanReingold), oo (OpenOrd)
	String layout
	Int? nsteps = 500
	Int? disk_space = 250
	Int? num_cpu = 2
	# memory in gigabytes
	Int? memory = 2

	call run_scrtools_graph_layout {
		input:
			input_file = input_file,
			output_name = output_name,
			layout = layout,
			nsteps = nsteps,
			disk_space = disk_space,
			num_cpu=num_cpu,
			memory=memory
	}
}

task run_scrtools_graph_layout {
    File input_file
    String output_name
    String layout
    Int? nsteps
    Int? disk_space
    Int? num_cpu
    String? memory

	command {
        set -e
        python <<CODE
        import h5py
        import pandas as pd
        import anndata
        import os
        import numpy as np
        from subprocess import check_call
        from scipy.sparse import csr_matrix
        # neighbors is at /uns/neighbors/connectivities
        # n_obs by n_obs
        f = h5py.File('${input_file}', 'r')
        n_obs = f['obs'].shape[0]
        group = f['uns/neighbors/connectivities']
        x = csr_matrix((group['data'][()], group['indices'][()], group['indptr'][()]), shape=(n_obs, n_obs))
        #obs_ids = f['obs']['index']
        writer = open('graph.net', 'w')
        writer.write('*Vertices ' + str(n_obs) + '\n')
        for i in range(n_obs):
            writer.write(str(i + 1) + ' "' + str(i + 1) + '"\n')
        writer.write('*Edges\n')
        rows, cols = x.nonzero()
        for i, j in zip(rows, cols):
            if j < i:
                writer.write(str(i + 1) + ' ' + str(j + 1) + ' ' + str(x[i, j]) + '\n')
        f.close()
        check_call(['java', '-Djava.awt.headless=true', '-Xmx${memory}g', '-cp', '/software/graph_layout/:/software/graph_layout/gephi-toolkit-0.9.2-all.jar', 'GraphLayout', 'graph.net', 'coords.txt', '${layout}', '${nsteps}'])

        df = pd.read_table('coords.txt', index_col=0)
        adata = anndata.read_h5ad('${input_file}', backed='r+')
        adata.obsm['X_draw_graph_' + '${layout}'] = np.array(list(zip(df.x, df.y)))
        adata.write('${output_name}')
        CODE
	}

	output {
		File output_h5ad = '${output_name}'
	}

	runtime {
		docker: "regevlab/scrtools_graphlayout"
		memory: "${memory} GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_cpu}"
		preemptible: 2
	}
}
