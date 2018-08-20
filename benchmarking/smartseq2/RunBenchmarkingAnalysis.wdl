import "BenchmarkingTasks.wdl" as analysis

workflow RunBenchmarkingAnalysis {
  
  meta {
    description: "Run SmartSeq2 Benchmarking pipeline. 4 modules are included in this pipeline.QC test, Comparative tests, Reproducibility test and Confounding factors test."
  }

  File base_datafile
  File updated_datafile
  File output_name
  File gtf_file
  File metadata_file
  File base_metrics
  File updated_metrics
  String metadata_keys
  String groups
  String low_cut
  String high_cut
  String met_keys
  Int    npcs
  String docker
  
  parameter_meta {
    base_datafile: "data matrix,count  or TPM matrix, of one pipeline"
    updated_datafile: "data matrix, count or TPM matrix, of second pipeline"
    base_metrics: "QC metrics, each row represents a metric and each column represents a cell"
    updated_metrics: "QC metrics, each row represents a metric and each column represents a cell"
    met_keys: "names of QC metrics to extract and compare"
    output_name: "output's prefix"
    low_cut: "low threshold cutoff to be used to filter out cells"
    high_cut: "high threshold cutoff to be used to filter out cells"
    metadata_file: "meta file"
    metadata_keys: "keys in metadata to be used to used as biological labels, such as cell type, cell lineage"
    groups: "labels to be used to identify conditions, such single cell vs bulk, donor1 vs donor2, control vs dieased"
    gtf_file: "gene annotaiton file, in gtf format"  
    npcs: "number of PCs to extract"
  } 

  call analysis.RunQCMetricsAnalysis as AnalysisQC {
    input:
      base_metrics = base_metrics,
      updated_metrics = updated_metrics,
      output_name = output_name,
      met_keys = met_keys,
      docker = docker
  }
  call analysis.RunComparativeAnalysis as ComparativeAnalysis {
    input:
      base_datafile = base_datafile,
      updated_datafile = updated_datafile,
      output_name = output_name,
      metadata_file = metadata_file,
      metadata_keys = metadata_keys,
      docker = docker
  }
  
 call analysis.RunGeneQuantificationAnalysis as CompareQuantification {
    input:
      base_datafile = base_datafile,
      updated_datafile = updated_datafile,
      output_name = output_name,
      gtf_file = gtf_file,
      low_cut = low_cut,
      high_cut = high_cut,
      docker = docker
   }
  
 call analysis.RunReproducibilityAnalysis  as AnalysisReproducibility {
    input: 
      base_datafile = base_datafile,
      updated_datafile = updated_datafile,
      output_name = output_name,
      gtf_file = gtf_file,
      metadata_file = metadata_file,
      groups = groups,
      docker = docker
  } 
 
  call analysis.RunConfoundingFactorAnalysis as AnalysisConfounding {
    input:
      base_datafile = base_datafile,
      updated_datafile = updated_datafile,
      base_metrics = base_metrics,
      updated_metrics = updated_metrics,
      output_name = output_name,
      npcs = npcs,
      metadata_file = metadata_file,
      meta_keys = metadata_keys,
      met_keys = met_keys,
      docker = docker
  }
  
 output {
  File qc_html = AnalysisQC.html
  File qc_test_stats = AnalysisQC.test_qc
  File confounding_html = AnalysisConfounding.html
  File confounding_results = AnalysisConfounding.test_confoundingfactors
  File reproducibility_html = AnalysisReproducibility.html
  File reproducibility_ttest_base = AnalysisReproducibility.ttest_base
  File reproducibility_ttest_updated = AnalysisReproducibility.ttest_updated
  File quantification_html = CompareQuantification.html
  File comparative_html = ComparativeAnalysis.html
  }
}
