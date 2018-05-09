import "BenchmarkingTasks.wdl" as analysis

workflow RunBenchmarkingAnalysis {
  meta {
    description: "Run SmartSeq2 Benchmarking pipeline. 4 modules are included in this pipeline.QC test, Comparative tests, Reproducibility test and Confounding factors test."
  }
  String base_datafile
  String updated_datafile
  String output_name
  String gtf_file
  String metadata_file
  String metadata_keys
  String met_keys
  String base_metrics
  String updated_metrics
  Int    npcs
  
  call analysis.RunQCMetricsAnalysis {
    input:
      base_metrics = base_metrics,
      updated_metrics = updated_metrics,
      output_name = output_name,
      met_keys = met_keys
  }

  call analysis.RunComparativeAnalysis {
    input:
      base_datafile = base_datafile,
      updated_datafile = updated_datafile,
      output_name = output_name,
      metadata_file = metadata_file,
      metadata_keys = metadata_keys
  }

  call analysis.RunReproducibilityAnalysis {
    input: 
      base_datafile = base_datafile,
      updated_datafile = updated_datafile,
      output_name = output_name,
      gtf_file = gtf_file,
      metadata_file = metadata_file
  } 
 
  call analysis.RunConfoundingFactorAnalysis {
    input:
      base_datafile = base_datafile,
      updated_datafile = updated_datafile,
      base_metrics = base_metrics,
      updated_metrics = updated_metrics,
      output_name = output_name,
      npcs = npcs,
      metadata_file = metadata_file,
      metaKeys = metadata_keys
  }
}
