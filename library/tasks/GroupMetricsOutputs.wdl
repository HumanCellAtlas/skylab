task GroupQCOutputs {
  Array[File] picard_row_outputs
  Array[File] picard_table_outputs
  File hisat2_stats
  File rsem_stats
  String output_name
  # Runtime
  String docker = "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.0.2"
  String mem = 2
  String cpu = 1
  String disks = 20 
  Int preemptible = 5
  Int max_retries = 0
  
  meta {
    description: "This task will group the Picard metrics"
  }
  parameter_meta {
    picard_outputs: "array of files generated by Picard"
    hisat2_stats: "statistics output of hisat2 alignment"
    rsem_stats: "statistics output of rsem "
    output_name: "name output files"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    max_retries: "(optional) retry this number of times if task fails -- use with caution, see skylab README for details"
  }
 command {
    set -e
    git clone --branch jx-ss2-bam-index https://github.com/HumanCellAtlas/skylab
    python skylab/pipelines/smartseq2_single_sample/GroupQCOutputs.py -f ${sep=' ' picard_row_outputs}  -t Picard -o Picard_group
    python skylab/pipelines/smartseq2_single_sample/GroupQCOutputs.py -f ${hisat2_stats} -t HISAT2 -o hisat2
    python skylab/pipelines/smartseq2_single_sample/GroupQCOutputs.py -f ${rsem_stats} -t RSEM -o rsem
    python skylab/pipelines/smartseq2_single_sample/GroupQCOutputs.py -f Picard_group.csv hisat2.csv rsem.csv -t Core -o "${output_name}_QCs"
    python skylab/pipelines/smartseq2_single_sample/GroupQCOutputs.py -f ${sep=' ' picard_table_outputs} -t PicardTable -o "${output_name}"
    ls *
}
  output{
    Array[File] group_files = glob("${output_name}_*.csv")
  }
  runtime {
    docker: docker
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 5
    maxRetries: 1
  }
}
