task checkOptimusInput {
  String chemistry;

  meta {
    description: "checks optimus input values and fails the pipeline immediately"
  }

  command {
    set -e
    
    ## Set pass to true
    pass=true

    ## Perform checks
    if [ ! ("${chemistry}" != "tenX_v2" && "${chemistry}" != "tenX_v3") ]
    then
	pass=1
	echo "ERROR: Invalid value \"${chemistry}\" for input \"chemistry\""
    fi

    ## fail if any tests failed
    if [ $pass ]
    then
      exit 0;
    else
      exit 1;
    fi
  }

  runtime {
    docker: "ubuntu:18.04"
    cpu: 1
    memory: "1.0 GB"
    disks: "local-disk 1 HDD"
  }
  
}