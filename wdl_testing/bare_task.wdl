
task BareTask {
  String user_input

  command <<<
    echo "user passed input: ${user_input}"
    echo "here are the environment variables $ENV"
  >>>
}
