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
    File output_a_file = write_lines([TaskA.output_string])
    String output_a_string = TaskA.output_string
  }
}

task TaskA {
  input {
    String string_input
  }

  command {
    echo ~{string_input}
  }

  output {
    String output_string = stdout()
  }

  runtime {
    docker: "alpine:latest"
    cpu: 1
    memory: "3.75 GiB"
    disk: "local-disk 10 HDD"
  }
}