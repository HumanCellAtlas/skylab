skylab_path="/Users/nbarkas/work/optimus/runs_Mar19/skylab"

cromwell-tools submit \
	       --username `vault_mint_cromwell_login` \
	       --password `vault_mint_cromwell_password` \
	       -w $skylab_path/pipelines/optimus/Optimus.wdl \
	       -d $skylab_path/library/tasks/* \
	       -l labels.json \
	       -o options.json \
	       --url https://cromwell.mint-dev.broadinstitute.org/ \
	       -i inputs_CB2_HiSeq_1.json |
               tee log.txt
