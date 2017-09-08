
task BareTask {
  String user_input

  command <<<
    echo "user passed input: ${user_input}"
    echo "here are the environment variables"
    TESTVAR="foo"
    echo $TESTVAR
  >>>

  output {
    Array[String] standard_out = read_lines(stdout())
  }
}
