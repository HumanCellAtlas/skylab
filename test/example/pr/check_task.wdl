task CheckPrintWorkflow {
    File output_file
    String expected_string

    command <<<

    observed_string=$(cat "${output_file}")
    if [ "${expected_string}" == "$observed_string" ]; then
        exit 0
    else
        exit 1
    fi

    >>>
}
