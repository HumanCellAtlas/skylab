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
    File output_a_file = read_lines(TaskC.output_string)
    String output_a_string = TaskC.output_string
  }
}

task TaskC {
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