# Running 10x Software in the Green Box

## Background

10x single cell data is one of the two pilot data types for the HCA Green Box, so
we need a workflow that processes 10x data to run on the Green Box infrastructure.

A reasonable initial workflow for 10x data analysis is the workflow released by
10x itself in a software suite called Cell Ranger. The Cell Ranger workflow hopefully
can serve as a decent starting point for iteration and improvement, and by using it
we get to import the experimentation and refinement that 10x has already done.

Cell Ranger is implemented in a way that emphasizes ease of use, so it presents a
very simple interface to users: a `cellranger` command that has a few subcommands.
The `cellranger` command launches a workflow that carries out some high-level
analysis task, such as producing a gene-by-cell matrix or converting BCLs to 
FASTQs. It runs this workflow locally or on an SGE-like cluster.

The Green Box can't use `cellranger` directly for two reasons. First, part of
the mission of the Green Box is workflow testing and improvement, so we want
finer control over components of the workflow, and we're willing to trade
away ease of use to get that. Second, the main execution environment for Green Box
workflows is Cromwell running on Google Cloud, which is unlike the execution
environments supported by `cellranger`.

So, we want to set aside the workflow execution and orchestration code that's run
by `cellranger` and work directly with the underlying analysis code. Unfortunately,
the underlying analysis code is implemented in an idiosycratic way that is tightly 
coupled with the orchestration code. To avoid having to rewrite or reimplement all
of the 10x analysis code for the pilot, we have to emulate a little bit of the
orchestration code so we can provide the interface that the underlying analysis
code is expecting.

## Martian

The 10x workflow engine is called Martian, and documentation for it is available
[here](http://martian-lang.org/). Martian has two entities: pipelines and stages.
Pipelines are composed of stages and generally accomplish large, complex tasks.
Stages are more focused and define more concrete tasks like "read alignment" or
"running bcl2fastq". Generally, `cellranger` only runs pipelines. But in the
Green Box, we want to be able to run stages directly, so we can do things like
replace an alignment tool or adjust its parameters.

Stages have three phases in Martian: "split", "main", and "join". These define
a scatter-gather pattern, where "split" breaks a task into parts, muliple "main"
instances perform the analysis, and "join" merges everything back together.
Some stages are simpler and only have a "main" phase.

Stages wrap a single python script that defines "split", "main", and "join"
functions. For example, the "EXTRACT_READS" stage wraps [this](https://github.com/10XGenomics/cellranger/blob/master/mro/stages/common/extract_reads/__init__.py)
python code. Looking at that script reveals some of the difficulty in running
the 10x analysis code directly. For example, the main function has two argument,
"args" and "outs" that we'd have to know how to construct. There are also a number
of global variables that the analysis code expects to be defined, and javascript-like
type system used in the stage definitions.

## martian_shell.py

Recently, 10x refactored Martian so it now includes a script called [martian_shell.py](https://github.com/martian-lang/martian/blob/master/adapters/python/martian_shell.py).
This can serve as a bridge between the bash command line environment that we have
when running WDL with Cromwwell and the environment expected by the 10x analysis
code.

`martian_shell.py` is executable and takes five arguments:
1. `stagecode_path` the path to the python code that defines the stage
2. `_run_type` split, main, or join
3. `metadata_path` not super relevant if we're only running one stage at a time
4. `files_path` path where some stages will write output files
5. `run_file` path where stages will write things like the current time

Also, each stage expects one or more json files to be present from which `martian_shell` will
read arguments.

### Split
Split phases expect a json file called "\_args" to be present that defines all the
members of the `args` object that the split function will need. So in the EXTRACT\_READS
example above, that's "chunks", "barcode_whitelist", and "initial_reads".

The split phase will write a file called "\_stage_defs" that is a json file with a key
"chunks" that points to a list of objects. Each object defines the scattered inputs
for subsequent "main" instances.

### Main
Main phases expect two json files to be present: "\_args" and "\_outs". The \_args file
is similar to the split phase; it needs all the members of `args` accessed in the main
function to be defined. Similarly, it requires the "\_outs" file to define all the members
of `outs` accessed in the main function.

`outs` are the outputs that the main phase is going to produce. So the \_outs file holds
either null values or future paths before the main phase is executed. In the EXTRACT\_READS
code above, the main phase will write files to `out.reads`, so in \_outs there should be
a `reads` key that defines the path to those future reads files. EXTRACT\_READS also
produces `outs.bam_comments` which is a list of strings. The value for this in \_outs
can just be null.

### Join
Join phases expect four files: \_args, \_outs, \_chunk\_defs, and \_chunk\_outs. The
first two are like above. \_chunk\_defs is the output of the split phase. \_chunk\_outs
is a list of the outputs of the main phases.

### \_jobinfo

Finally, there's a json file expected by all phases that's just hard coded. It's called
"\_jobinfo" and consists of this:

```
{
  "profile_mode": "disable",
  "stackvars_flag": "disable",
  "monitor_flag": "disable",
  "invocation": "whatever",
  "version": {
    "martian": "'HCA-v2.3.0-rc2'",
    "pipelines": "HCA-7fe84e2"
  }
}
```

It doesn't really matter, but you will see it in the WDL.

## WDLizing the Martian Interface

Martian uses a handful of json files to define inputs and outputs. But we don't want all
of our WDL task to have inputs like `File args`, since that would make the true interface
of the step completely opaque. And since the Cell Ranger code is just a starting point,
we want to know what the inputs and outputs of each step are so we can start to experiment
with and improve them.

To achieve this, we avoid passing around the Martian json files wherever possible in favor
of defining explicit inputs and outputs for each WDL task. The tradeoff is that within each
task, the Martian json files need to be constructed and parsed before and after running
the code for the stage.

Currently, this is all done with `jq`, a command line program with a concise syntax for
manipulating json files.

## Files in Martian

Martian assumes that there's a persistent filesystem shared among all the stages and phases,
so it isn't careful about defining file inputs and outputs. So a stage may create a file but
then appear to only output a string that's the path to the new file.

In an SGE environment with a shared filesystem, that would work fine since that string would
be a valid path to the file in any other stage. But in a cloud environment, that won't be
true, so we have to handle file staging and transfers on our own.

Also, when Cromwell has a file input, the path to that input is specific to the task execution.
It contains references to the task name, UUIDs, etc. Martian will sometimes try to store and
pass on that path to downstream stages where the path will be invalid.

To get around this, we often have to hard link input files to the working directory of the
task execution. This way, we can create paths that are valid across different tasks and
not have to worry about Martian trying to refer to some path that was only valid in a
different stage.
