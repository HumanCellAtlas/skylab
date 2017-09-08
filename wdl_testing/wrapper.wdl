import "bare_task.wdl" as bare_task

workflow WrapperForTask {
  String user_input

  call bare_task.BareTask {
    input:
     user_input = user_input
  }
}
