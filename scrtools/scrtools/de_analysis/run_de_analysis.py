from datetime import datetime
import scanpy.api as sc

from . import de_analysis

def run_de_analysis(input_file, output_name, threshold):
	adata = sc.read(input_file)

	print("Begin t_test.")
	print(datetime.now())
	de_analysis.t_test(adata)
	print(datetime.now())
	print("t_test is done.")

	excel_file = output_name + "_de_analysis_t.xlsx"
	de_analysis.write_results_to_excel(excel_file, adata, "t", threshold = threshold)
	print(excel_file + " is written.")

	print("Begin fisher_test.")
	print(datetime.now())
	de_analysis.fisher_test(adata)
	print(datetime.now())
	print("fisher_test is done.")

	excel_file = output_name + "_de_analysis_fisher.xlsx"
	de_analysis.write_results_to_excel(excel_file, adata, "fisher", threshold = threshold)
	print(excel_file + " is written.")

	adata.write(output_name + "_de.h5ad")
	print("Results are written.")
