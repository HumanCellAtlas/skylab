import "test_task.wdl"
import "check_task.wdl"

workflow Test {
    String printed_string

    call test_task.Print as print {
        input:
            printed_string = printed_string
    }

    call check_task.CheckPrintWorkflow {
        input:
            output_file = print.output_file,
            expected_string = printed_string
    }
}
