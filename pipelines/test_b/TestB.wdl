version 1.0

workflow TestB {

  input {
    String input_b
  }

  call TaskB {
    input:
      string_input = input_b
  }

  output {
    File output_a_file = TaskB.output_file
    String output_a_string = TaskB.output_string
  }
}

task TaskB {
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