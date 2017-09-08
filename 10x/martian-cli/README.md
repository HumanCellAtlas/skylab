## martian_cli

Alternative interface for the 10x Genomics Martian workflow executor.

10x provides a suite of tools for analysis of their single-cell data called Cell Ranger.
Cell Ranger contains the actual analysis code as well as a workflow executor called
Martian.

For the HCA Green Box, we want to be able to use the 10x analysis code, at least as a
starting point. But, we cannot use the workflow executor as it doesn't easily allow
modifications to the pipeline, and it's currently closed source.

`martian_cli` is an alternative way to run 10x analysis code that doesn't use Martian
and ideally integrates easily with community-driver workflow languages like WDL
and CWL.
