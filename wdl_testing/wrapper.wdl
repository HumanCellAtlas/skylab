import "bare_task.wdl" as bare_task

workflow WrapperForTask {
  String user_input

  call bare_task.BareTask as bt {
    input:
     user_input = user_input
  }

  output {
    Array[String] standard_output = bt.standard_out
  }
}
