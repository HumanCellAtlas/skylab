version 1.0

workflow SnapATACSecondaryAnalysis {
    input {
        File snap_file
        String sample_name
    }

    parameter_meta {
        snap_file: "snap file"
        sample_name: "sample name"
    }

    call RunAnalysis {
        input:
            snap_file = snap_file,
            sample_name = sample_name
    }

    output {
        # fixme: use array of files
        File output_barcode_pdf = RunAnalysis.output_barcode_pdf
        # File output_bin_size_png = RunAnalysis.output_bin_size_png
        File output_bin_coverage_pdf = RunAnalysis.output_bin_coverage_pdf
        File output_pca_elbow_pdf = RunAnalysis.output_pca_elbow_pdf
        File output_pca_pw_pdf = RunAnalysis.output_pca_pw_pdf
        File output_tsne_pdf = RunAnalysis.output_tsne_pdf
        File output_umap_pdf = RunAnalysis.output_umap_pdf
        File output_marker_genes_pdf = RunAnalysis.output_marker_genes_pdf
        Array[File] output_macs2 = RunAnalysis.output_macs2
    }
}

task RunAnalysis {
    input {
        File snap_file
        String sample_name
        String path_outputs = "/opt/outputs"
        String docker_image = "hisplan/snap-atac-secondary-analysis:0.0.7-alpha"
    }

    Int num_threads = 16

    command <<<
        set -euo pipefail

        mkdir -p ~{path_outputs}
        # Rscript /opt/secondary-analysis.R -o ~{path_outputs} -i ~{snap_file} -d /opt/data -s ~{sample_name}
        Rscript /opt/secondary-analysis.R ~{snap_file} ~{path_outputs} /opt/data ~{sample_name}

    >>>

    output {
        # fixme: use an array of files with glob
        File output_barcode_pdf = "barcode.pdf"
        # File output_bin_size_png = "bin-size.png"
        File output_bin_coverage_pdf = "bin-coverage.pdf"
        File output_pca_elbow_pdf = "pca-elbow.pdf"
        File output_pca_pw_pdf = "pca-pw.pdf"
        File output_tsne_pdf = "tsne.pdf"
        File output_umap_pdf = "umap.pdf"
        File output_marker_genes_pdf = "marker-genes.pdf"
        Array[File] output_macs2 = glob("macs2_out.*")
    }

    runtime {
        docker: docker_image
        # fixme:
        disks: "local-disk 500 HDD"
        cpu: num_threads
        memory: "32 GB"
    }
}
