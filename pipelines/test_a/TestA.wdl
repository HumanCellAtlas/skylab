version 1.0

workflow TestA {

  input {
    String input_a
  }

  call TaskA {
    input:
      string_input = input_a
  }

  output {
    File output_a_file = TaskA.output_file
    String output_a_string = TaskA.output_string
  }
}

task TaskA {
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