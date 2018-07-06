## Releasing Approved HCA Pipelines

### AWG Approval Process

...

### Publishing in methods repositories
In order to make our pipelines more widely available, we add them to methods repositories so potential users can discover and execute them. Currently, we publish our pipelines in [Dockstore](https://dockstore.org/), which implements the GA4GH [Tool Registry Service](https://github.com/ga4gh/tool-registry-schemas) API.

#### Registering a workflow in Dockstore
Dockstore has some expectations about the organization of github repos that do not quite match the way `skylab` is organized. It also expects input files to be http urls, not gs urls. So, there are a couple things to do within skylab before the workflow can be registered in Dockstore:

<dl>
  <dt>Switch to the <code>dockstore</code> branch and rebase it on <code>master</code> (or the release branch)</dt>
  <dd>We want to keep Dockstore-specific detail outside of <code>master</code> so things don't look confusing for people looking at the repo. And, Dockstore allows workflows to be published from specific branches of a repo, so we put everything specific to Dockstore in a branch that we keep up-to-date with released pipeline code.</dd>
  <dt>Add symlinks to imported WDLs adjacent to the entry point WDL</dt>
  <dd>We have WDLs that import other WDLs (that import other WDLs...), and we rely on our submission code to handle finding those dependencies and staging them to the imports all resolve correctly. Dockstore has no way of knowing about this though, so we need to insert links into the repo so Dockstore can find the dependencies.
  Our imports never specify any paths, so we expect all dependency WDLs to be in a flat directory along with the main, entry point WDL. So, we just need to add symlinks to all the imported WDLs alongside the WDLs we're actually going to run. For example:
  <pre>
  pipelines/smartseq2_single_sample/SmartSeq2SingleSample.wdl
  pipelines/smartseq2_single_sample/HISAT.wdl (link to ../../library/tasks/HISAT.wdl)
  pipelines/smartseq2_single_sample/Picard.wdl (link to ../../library/tasks/Picard.wdl)
  ...</pre></dd>
  <dt>Prepare Dockstore test input files</dt>
  <dd>Our current test input files specify gs urls. One day these will be HCA DSS urls or maybe even GA4GH DOS urls, but that day is not today. Dockstore wants http urls, so we need to prepare a new input file:
  <pre>
  {
    "SmartSeq2SingleCell.gtf_file": "http://storage.googleapis.com/hca-dcp-mint-test-data/reference/GRCh38_Gencode/gencode.v27.primary_assembly.annotation.gtf",
    "SmartSeq2SingleCell.hisat2_ref_trans_name": "gencode_v27_trans_rsem",
    "SmartSeq2SingleCell.rrna_intervals": "http://storage.googleapis.com/hca-dcp-mint-test-data/reference/GRCh38_Gencode/gencode.v27.rRNA.interval_list",
    ...</pre></dd>
</dl>

Now we're ready to go to [Dockstore](https://dockstore.org/) and actually register the pipeline:

- Go to dockstore.org and login in via github. You'll want to use the github account you use to commit to `skylab`
- Click on "My Workflows", then "Register Workflow", then enter details about the entry point WDL. The repo name should be `HumanCellAtlas/skylab`, and the workflow path should be something like `/pipelines/smartseq2_single_sample/ss2_single_sample.wdl`.
- Click "Register Workflow". If nothing happens, wait a few seconds and click "Register Workflow" again. This should trigger an error, and then you can click "Close"
- Click on "Versions" and select the `dockstore` branch radio button.
- Click "Publish". When you click on "Files", you should see the entry point WDL and all of the dependency WDLs listed. If you don't see them, it's possible that the selected branch has reverted to the most recently updated branch rather than `dockstore`. So make sure `dockstore` is selected from the branch dropdown menu.

TODO: Verify that users can run the workflow. Presently, none of the suggested commands or launch buttons will work. That should change in the future.
