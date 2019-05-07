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
    File output_a_file = read_lines(TaskB.output_string)
    String output_a_string = TaskB.output_string
  }
}

task TaskB {
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