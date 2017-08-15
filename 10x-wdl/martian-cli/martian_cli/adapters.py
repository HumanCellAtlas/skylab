import json
import os
import sys

import martian

from martian_cli import common


def _construct_outs(martian_io_fields):
    """Create an outs object for Martian functions.

    Main and join phases have an input parameter called "outs". This seems to be a record with
    keys for the expected outputs of the phase and values of None if the output is not a file
    or a string that's the same as the key + the file type extension if it is a file.
    """
    outs = {}
    for field in martian_io_fields:
        if field.type in common.MARTIAN_FILETYPES:
            outs[field.name] = field.name + '.' + field.type
        elif field.type == "path":
            outs[field.name] = field.name
        else:
            outs[field.name] = None
    return martian.Record(outs)

def _write_outs(outs, martian_io_fields):
    """Write outputs of a phase to files.

    Martian functions have a return value, but those get lost when they're called as
    a WDL stage. So the Martian CLI needs to write outputs to files that the WDL runner
    can read.
    """

    for field in martian_io_fields:
        # Maps require special handling so we can use the WDL read_map function
        if field.type == 'map':
            with open(field.name, 'w') as outf:
                for key, value in getattr(outs, field.name).items():
                    outf.write("\t".join([str(key), str(value)]))
                    outf.write("\n")
        elif field.type not in common.MARTIAN_FILETYPES and field.type != "path":
            with open(field.name, 'w') as outf:
                json.dump(getattr(outs, field.name), outf)

def _common(stage, phase):
    """Martian functions expect a bunch of global scope variables that behave in a wildly
    unpythonic way. That can be mimicked by shamefully sticking things into __builtin__.
    """
    import __builtin__
    __builtin__.metadata = martian.Metadata(stage.source, "files", "run", phase)
    sys.path.append(os.path.dirname(stage.source))
    __builtin__.module = __import__(os.path.basename(stage.source))

def split(stage, cli_args):
    """Run a split phase of a stage."""

    _common(stage, "split")

    args = martian.Record(
        {k.name: getattr(cli_args, k.name) for k in stage.inputs})

    # Now continue with the theme of screwing with the environment, and exec the module
    env = globals()
    env.update(locals())
    exec("stage_defs = module.split(args)", env, env)

    _write_outs(stage_defs, stage.splits)

    return stage_defs["chunks"]

def main(stage, cli_args):
    """Run a main phase of a stage."""

    _common(stage, "main")

    # Build the args. The args of a main stage are all the "inputs" plus everything
    # in "split using"
    args = martian.Record(
        {k.name: getattr(cli_args, k.name) for k in stage.inputs + stage.splits})

    # Build the outs object, which is just the stage output fields
    outs = _construct_outs(stage.outputs)

    # Now continue with the theme of screwing with the environment, and exec the module
    env = globals()
    env.update(locals())
    exec("module.main(args, outs)", env, env)
    print "outs", outs
    _write_outs(outs, stage.outputs)

    return outs

def join(stage, cli_args):
    """Run a join phase of a stage."""

    _common(stage, "join")

    # The inputs to a join stage are complicated. There are four inputs
    # 1. args - these are the same as the args to split, which is just stage inputs
    # 2. outs - this is the stage outputs, constructed in the usual way
    # 3. chunk_defs - a list of the inputs provided to each of the main steps. So these
    #        are fields from splits
    # 4. chunk_outs - a list of the outputs from each of the main steps. So these are
    #        fields from outs

    args = martian.Record(
        {k.name: getattr(cli_args, k.name) for k in stage.inputs})
    outs = _construct_outs(stage.outputs)

    all_splits = {k.name: getattr(cli_args, k.name + "_split") for k in stage.splits}
    chunk_defs = [martian.Record(dict(zip(all_splits.keys(), v))) for v in zip(*(all_splits.values()))]
    all_chunks = {k.name: getattr(cli_args, k.name + "_output") for k in stage.outputs}
    chunk_outs = [martian.Record(dict(zip(all_chunks.keys(), v))) for v in zip(*(all_chunks.values()))]

    print "args", args
    print "outs", outs
    print "chunk_defs", chunk_defs
    print "chunk_outs", chunk_outs

    # Now continue with the theme of screwing with the environment, and exec the module
    env = globals()
    env.update(locals())
    exec("module.join(args, outs, chunk_defs, chunk_outs)", env, env)

    _write_outs(outs, stage.outputs)

    return outs
