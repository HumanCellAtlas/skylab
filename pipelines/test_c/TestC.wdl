version 1.0

workflow TestC {

  input {
    String input_c
  }

  call TaskC {
    input:
      string_input = input_c
  }

  output {
    File output_a_file = TaskC.output_file
    String output_a_string = TaskC.output_string
  }
}

task TaskC {
  input {
    String string_input
  }

  command {
    echo ~{string_input}
    echo ~{string_input} > output_file
  }

  output {
    String output_string = read_lines(stdout())[0]
    File output_file = "output_file"
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
    cpu: 1
    memory: "3.75 GiB"
    disk: "local-disk 10 HDD"
  }
}