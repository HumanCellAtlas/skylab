# Task Description
# ================
#
# This task simply prints an input string "printed_string" to file and exposes that
# file as "output_file".

task Print {

    String printed_string

    command <<<

        # print the string to file
        echo "${printed_string}" >> output.txt

    >>>

    # need not request disks for this; start-up disk is adequate
    runtime {
        docker: "ubuntu:bionic-20190307"
        cpu: 1
        memory: "3.75 GB"
    }

    output {
        File output_file = "output.txt"
    }
}
