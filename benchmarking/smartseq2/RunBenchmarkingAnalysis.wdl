import "run_combine_QCs.wdl" as CombineQC
import "run_matrix_comparison.wdl" as CompareDataMatrix
import "run_matrix_reproducibility.wdl" as ReproduciblityDataMatrix
import "run_qc_metrics_analysis.wdl" as CompareQC
import "run_confounding_factor_analysis.wdl" as ConfoundingDataMatrix

workflow RunBenchmarkingAnalysis {
  String scripts_dir
  String base_datafile
  String updated_datafile
  String output_name
  String gtf_file
  String metadata_file
  String metKeys
  String base_metrics
  String updated_metrics
  Int    npcs
  Int    nthreshold

  call CombineQC.RunCombineQC as CombineBase {
    input:
      scripts_dir = scripts_dir,
      datafile = base_datafile,
      metfile = base_metrics,
      output_name = output_name,
      gtf_file = gtf_file,
      nthreshold = nthreshold
  }

  call CombineQC.RunCombineQC as CombineUpdates {
    input:
      scripts_dir = scripts_dir,
      datafile = updated_datafile,
      metfile = updated_metrics,
      output_name = output_name,
      gtf_file = gtf_file,
      nthreshold = nthreshold
  }
  
  call CompareQC.RunQCAnalysis {
    input:
      scripts_dir = scripts_dir,
      base_qc = CombineBase.CombinedQC,
      updated_qc = CombineUpdates.CombinedQC,
      output_dir = output_name,
      metrics_keys = metKeys
  }

  call CompareDataMatrix.RunDataMatrixComparison {
    input:
      scripts_dir = scripts_dir,
      base_datafile = base_datafile,
      updated_datafile = updated_datafile,
      output_name = output_name,
      metadata_file = metadata_file
  }

  call ReproduciblityDataMatrix.RunDataMatrixReproducibility {
    input: 
      base_datafile = base_datafile,
      updated_datafile = updated_datafile,
      output_name = output_name,
      gtf_file = gtf_file,
      scripts_dir = scripts_dir,
      metadata_file = metadata_file
  } 
 
  call ConfoundingDataMatrix.RunConfoundingFactorsAnalysis {
    input:
      base_datafile = base_datafile,
      updated_datafile = updated_datafile,
      base_metrics = CombineBase.CombinedQC,
      updated_metrics = CombineUpdates.CombinedQC,
      output_name = output_name,
      scripts_dir = scripts_dir,
      npcs = npcs
  }
}
