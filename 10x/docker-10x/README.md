## 10x Docker Image

Running 10x WDLs with Cromwell is best done within a docker container that provides
an environment with Cell Ranger and its dependencies. The image defined here does
that. It does three basic things:

1. Install bcl2fastq from Illumina. Cell Ranger requires this to be in the path when
running the "make fastq" pipeline.
2. Extract the Cell Ranger tarball. This is available from 10x and contains all of the
10x code that the WDL will want to run.
3. Install the martian\_cli wrapper. This allows the WDL to call 10x analysis code
directly without needing 10x's Martian workflow executor.
4. Add a little "run_in_10x" script that sources the Cell Ranger environment before
running a command.

A couple notes:
1. Only Cell Ranger 2.0.0 works with martian\_cli. Cell Ranger 2.0.1 definitely does
not work.
2. Normally the "run_in_10x" script would be an entrypoint, but Cromwell does not
support that.
