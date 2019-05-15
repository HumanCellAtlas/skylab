# ENCODE DCC RNA-seq pipeline
# https://github.com/ENCODE-DCC/rna-seq-pipeline/tree/v1.0

workflow rna {
    meta {
        description: "Analyze bulk RNA-seq data using the ENCODE DCC RNA-seq pipeline."
    }
    # endedness: paired or single
    String endedness
    # fastqs_R1: fastq.gz files for Read1 (only these if single-ended)
    Array[File] fastqs_R1
    # fastqs_R2: fastq.gz files for Read2 (omit if single-ended) in order
    # corresponding to fastqs_R1
    Array[File] fastqs_R2 = []
    # aligner: star for now, more added if/when needed
    String aligner
    # bamroot: root name for output bams. For example foo_bar will
    # create foo_bar_genome.bam and foo_bar_anno.bam
    String bamroot
    # strandedness: is the library strand specific (stranded or unstranded)
    String strandedness
    # strandedness_direction (forward, reverse, unstranded)
    String strandedness_direction
    # chrom_sizes: chromosome sizes file
    File chrom_sizes

    ## task level variables that are defined globally to make them visible to DNANexus UI

    # ALIGN
    # index: aligner index (tar.gz)
    File align_index
    Int align_ncpus
    Int align_ramGB
    # indexdir: where to extract the star index, relative to cwd
    String? indexdir
    # libraryid: identifier which will be added to bam headers
    String? libraryid
    String? align_disk

    # KALLISTO

    Int kallisto_number_of_threads
    Int kallisto_ramGB
    File kallisto_index
    Int? kallisto_fragment_length
    Float? kallisto_sd_of_fragment_length
    String? kallisto_disk

    # BAM_TO_SIGNALS

    Int bam_to_signals_ncpus
    Int bam_to_signals_ramGB
    String? bam_to_signals_disk

    # RSEM_QUANT

    # rsem_index: location of the RSEM index archive (tar.gz)
    File rsem_index
    # rnd_seed: random seed used for rsem
    Int rnd_seed = 12345
    Int rsem_ncpus
    Int rsem_ramGB
    String? rsem_disk

    # RNA_QC

    File rna_qc_tr_id_to_gene_type_tsv
    String? rna_qc_disk

    # MAD_QC

    String? mad_qc_disk

    # RUNTIME

    Int? preemptible

    parameter_meta {
        sample_id: "Name of sample to run CellRanger count on"
        endedness: "paired or single"
        fastqs_R1: "fastq.gz files for Read1 (only these if single-ended)"
        fastqs_R2: "fastq.gz files for Read2 (omit if single-ended) in order corresponding to fastqs_R1"
        aligner: "star for now, more added if/when needed"
        bamroot: "root name for output bams. For example foo_bar will create foo_bar_genome.bam and foo_bar_anno.bam"
        strandedness: "is the library strand specific (stranded or unstranded)"
        strandedness_direction: "forward, reverse, unstranded"
        chrom_sizes: "chromosome sizes file"
        indexdir: "where to extract the star index, relative to cwd"
        libraryid: "identifier which will be added to bam headers"
        rsem_index: "location of the RSEM index archinve (tar.gz)"
        rnd_seed: "random seed used for rsem"
        preemptible: "optional number of preemptible tries all tasks will take"
    }

    Array[Array[File]] fastqs_ = if length(fastqs_R2)>0 then transpose([fastqs_R1, fastqs_R2]) else transpose([fastqs_R1])

    scatter (i in range(length(fastqs_))) {
        call align { input:
            endedness = endedness,
            fastqs = fastqs_[i],
            index = align_index,
            aligner = aligner,
            indexdir = indexdir,
            libraryid = libraryid,
            bamroot = "rep"+(i+1)+bamroot,
            ncpus = align_ncpus,
            ramGB = align_ramGB,
            disks = align_disk,
            preemptible = preemptible
        }

        call bam_to_signals { input:
            input_bam = align.genomebam,
            chrom_sizes = chrom_sizes,
            strandedness = strandedness,
            bamroot = "rep"+(i+1)+bamroot+"_genome",
            ncpus = bam_to_signals_ncpus,
            ramGB = bam_to_signals_ramGB,
            disks = bam_to_signals_disk,
            preemptible = preemptible
        }

        call rsem_quant { input:
            rsem_index = rsem_index,
            rnd_seed = rnd_seed,
            anno_bam = align.annobam,
            endedness = endedness,
            read_strand = strandedness_direction,
            ncpus = rsem_ncpus,
            ramGB = rsem_ramGB,
            disks = rsem_disk,
            preemptible = preemptible
        }
    }

    scatter (i in range(length(fastqs_))) {
        call kallisto { input:
            fastqs = fastqs_[i],
            endedness = endedness,
            strandedness_direction = strandedness_direction,
            kallisto_index = kallisto_index,
            number_of_threads = kallisto_number_of_threads,
            ramGB = kallisto_ramGB,
            fragment_length = kallisto_fragment_length,
            sd_of_fragment_length = kallisto_sd_of_fragment_length,
            disks = kallisto_disk,
            out_prefix = "rep"+(i+1)+bamroot,
            preemptible = preemptible
        }
    }

    # if there are exactly two replicates, calculate the madQC metrics and draw a plot

    if (length(fastqs_R1) == 2) {
        call mad_qc { input:
            quants1 = rsem_quant.genes_results[0],
            quants2 = rsem_quant.genes_results[1],
            disks = mad_qc_disk,
            preemptible = preemptible
        }
    }

    scatter (i in range(length(align.annobam))) {
        call rna_qc { input:
            input_bam = align.annobam[i],
            tr_id_to_gene_type_tsv = rna_qc_tr_id_to_gene_type_tsv,
            output_filename = "rep"+(i+1)+bamroot+"_qc.json",
            disks = rna_qc_disk,
            preemptible = preemptible
        }
    }

    output {
      Array[File] anno_flagstat = align.anno_flagstat
      Array[File] annobam = align.annobam
      Array[File] genome_flagstat = align.genome_flagstat
      Array[File] genomebam = align.genomebam
      Array[File] align_log = align.log
      Array[File] align_python_log = align.python_log
      Array[Array[File]] bam_to_signals_all = bam_to_signals.all
      Array[File] bam_to_signals_python_log = bam_to_signals.python_log
      Array[Array[File]] bam_to_signals_unique = bam_to_signals.unique
      Array[File] kallisto_python_log = kallisto.python_log
      Array[File] kallisto_quants = kallisto.quants
      File? madQCmetrics = mad_qc.madQCmetrics
      File? madQCplot = mad_qc.madQCplot
      File? mad_qc_python_log = mad_qc.python_log
      Array[File] rna_qc_python_log = rna_qc.python_log
      Array[File] rnaQC = rna_qc.rnaQC
      Array[File] genes_results = rsem_quant.genes_results
      Array[File] isoforms_results = rsem_quant.isoforms_results
      Array[File] number_of_genes = rsem_quant.number_of_genes
      Array[File] rsem_quant_python_log = rsem_quant.python_log
    }
}


    ## tasks
    task align {
        Array[File] fastqs
        String endedness
        String aligner
        File index
        String? indexdir
        String? libraryid
        String bamroot
        Int ncpus
        Int ramGB
        String? disks
        Int? preemptible

        command {
            python3 $(which align.py) \
                ${if length(fastqs)<2 then "--fastqs " + fastqs[0] else "--fastqs " + fastqs[0] + " " + fastqs[1]} \
                --endedness ${endedness} \
                --aligner ${aligner} \
                --index ${index} \
                ${"--indexdir " + indexdir} \
                ${"--libraryid " + libraryid} \
                ${"--bamroot " + bamroot} \
                ${"--ncpus " + ncpus} \
                ${"--ramGB " + ramGB}
        }

        output {
            File genomebam = glob("*_genome.bam")[0]
            File annobam = glob("*_anno.bam")[0]
            File genome_flagstat = glob("*_genome_flagstat.txt")[0]
            File anno_flagstat = glob("*_anno_flagstat.txt")[0]
            File log = glob("*_Log.final.out")[0]
            File python_log = glob("align.log")[0]
        }

        runtime {
          cpu: ncpus
          memory: "${ramGB} GB"
          disks : select_first([disks, "local-disk 100 SSD"])
          docker : "quay.io/encode-dcc/rna-seq-pipeline:v1.0"
          preemptible: select_first([preemptible, 5])
          bootDiskSizeGb: 10
          noAddress: false
        }
    }

    task bam_to_signals {
        File input_bam
        File chrom_sizes
        String strandedness
        String bamroot
        Int ncpus
        Int ramGB
        String? disks
        Int? preemptible


        command {
            python3 $(which bam_to_signals.py) \
                --bamfile ${input_bam} \
                --chrom_sizes ${chrom_sizes} \
                --strandedness ${strandedness} \
                --bamroot ${bamroot}
        }

        output {
            Array[File] unique = glob("*niq.bw")
            Array[File] all = glob("*ll.bw")
            File python_log = glob("bam_to_signals.log")[0]
        }

        runtime {
            cpu: ncpus
            memory: "${ramGB} GB"
            disks : select_first([disks,"local-disk 100 SSD"])
            docker : "quay.io/encode-dcc/rna-seq-pipeline:v1.0"
            preemptible: select_first([preemptible, 5])
            bootDiskSizeGb: 10
            noAddress: false
        }
    }

    task rsem_quant {
        File rsem_index
        File anno_bam
        String endedness
        String read_strand
        Int rnd_seed
        Int ncpus
        Int ramGB
        String? disks
        Int? preemptible

        command {
            python3 $(which rsem_quant.py) \
                --rsem_index ${rsem_index} \
                --anno_bam ${anno_bam} \
                --endedness ${endedness} \
                --read_strand ${read_strand} \
                --rnd_seed ${rnd_seed} \
                --ncpus ${ncpus} \
                --ramGB ${ramGB}
        }

        output {
            File genes_results = glob("*.genes.results")[0]
            File isoforms_results = glob("*.isoforms.results")[0]
            File python_log = glob("rsem_quant.log")[0]
            File number_of_genes = glob("*_number_of_genes_detected.json")[0]
        }

        runtime {
            cpu: ncpus
            memory: "${ramGB} GB"
            disks : select_first([disks,"local-disk 100 SSD"])
            docker : "quay.io/encode-dcc/rna-seq-pipeline:v1.0"
            preemptible: select_first([preemptible, 5])
            bootDiskSizeGb: 10
            noAddress: false
        }
    }

    task kallisto {
        Array[File] fastqs
        File kallisto_index
        String endedness
        String strandedness_direction
        Int number_of_threads
        Int ramGB
        String out_prefix
        Int? fragment_length
        Float? sd_of_fragment_length
        String? disks
        Int? preemptible

        command {
            python3 $(which kallisto_quant.py) \
                --fastqs ${sep=' ' fastqs} \
                --number_of_threads ${number_of_threads} \
                --strandedness ${strandedness_direction} \
                --path_to_index ${kallisto_index} \
                --endedness ${endedness} \
                ${"--fragment_length " + fragment_length} \
                ${"--sd_of_fragment_length " + sd_of_fragment_length} \
                ${"--out_prefix " + out_prefix}
        }

        output {
            File quants = glob("kallisto_out/*_abundance.tsv")[0]
            File python_log = glob("kallisto_quant.log")[0]
        }

        runtime {
            cpu: number_of_threads
            memory: "${ramGB} GB"
            disks: select_first([disks, "local-disk 100 SSD"])
            docker : "quay.io/encode-dcc/rna-seq-pipeline:v1.0"
            preemptible: select_first([preemptible, 5])
            bootDiskSizeGb: 10
            noAddress: false
        }
    }

    task mad_qc {
    File quants1
    File quants2
    String? disks
    Int? preemptible

        command {
            python3 $(which mad_qc.py) \
                --quants1 ${quants1} \
                --quants2 ${quants2} \
                --MAD_R_path $(which MAD.R)
        }

        output {
            File madQCplot = glob("*_mad_plot.png")[0]
            File madQCmetrics = glob("*_mad_qc_metrics.json")[0]
            File python_log = glob("mad_qc.log")[0]
        }

        runtime {
            cpu: 2
            memory: "3400 MB"
            disks: select_first([disks,"local-disk 100 SSD"])
            docker : "quay.io/encode-dcc/rna-seq-pipeline:v1.0"
            preemptible: select_first([preemptible, 5])
            bootDiskSizeGb: 10
            noAddress: false
        }
    }

    task rna_qc {
        File input_bam
        File tr_id_to_gene_type_tsv
        String output_filename
        String? disks
        Int? preemptible

        command {
            python3 $(which rna_qc.py) \
                --input_bam ${input_bam} \
                --tr_id_to_gene_type_tsv ${tr_id_to_gene_type_tsv} \
                --output_filename ${output_filename}
        }

        output {
            File rnaQC = glob("*_qc.json")[0]
            File python_log = glob("rna_qc.log")[0]
        }

        runtime {
            cpu: 2
            memory: "1024 MB"
            disks: select_first([disks, "local-disk 100 SSD"])
            docker : "quay.io/encode-dcc/rna-seq-pipeline:v1.0"
            preemptible: select_first([preemptible, 5])
            bootDiskSizeGb: 10
            noAddress: false
        }
    }