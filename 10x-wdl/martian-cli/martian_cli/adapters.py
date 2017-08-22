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

def _write_chunks(stage_defs, chunks_object_filename):
    """Write the chunk outputs of a split step to json files.

    Martian split steps tend to specify multiple fields that you're supposed to scatter
    over. But, WDL only lets you scatter over one array, or two if you use zip. So we
    just have to shove everything into json files and scatter over those.
    """
    chunks = stage_defs["chunks"]

    for idx, chunk in enumerate(chunks):
        with open("martian_split_{}".format(idx), 'w') as split_file:
            json.dump(chunk, split_file)

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

    _write_chunks(stage_defs, "martian_splits")

    return stage_defs

def main(stage, cli_args):
    """Run a main phase of a stage."""

    _common(stage, "main")

    # If there's a "split_file" argument, that's the json file that was created
    # by the split phase. We need to try to parse that and insert the values
    # into the args Record that the main function is going to get.
    cli_dict = vars(cli_args)
    if hasattr(cli_args, "split_file") and cli_args.split_file:
        split_dict = json.load(open(cli_args.split_file))
        cli_dict.update(split_dict)

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

    # Args and outs are simple
    args = martian.Record(
        {k.name: getattr(cli_args, k.name) for k in stage.inputs})
    outs = _construct_outs(stage.outputs)

    # Now handle the chunk_defs and chunk_outs

    # Read all the arguments that are passed from the split step.
    all_splits = {k.name: getattr(cli_args, k.name + "_split") for k in stage.splits}
    # Filter out any omitted split values.
    all_splits = {k: v for k, v in all_splits.items() if v is not None}
    # Transpose into a list of dicts rather than dict of lists, and create martian.Records
    chunk_defs = [martian.Record(dict(zip(all_splits.keys(), v))) for v in zip(*(all_splits.values()))]

    # Now similarly, read arguments that are outputs from main steps
    all_chunks = {k.name: getattr(cli_args, k.name + "_output") for k in stage.outputs}
    all_chunks = {k: v for k, v in all_chunks.items() if v is not None}

    # Now transpose but don't put into martian.Records yet
    flattened_chunks = [dict(zip(all_chunks.keys(), v)) for v in zip(*(all_chunks.values()))]

    # The outputs from the main steps might be json files. If they are, then parse them and insert the
    # values into the dicts in place of the file name
    parsed_chunks = []
    for flat_chunk in flattened_chunks:
        parsed_chunk = {}
        for key in flat_chunk:
            try:
                parsed_json = json.load(open(flat_chunk[key]))
                parsed_chunk[key] = parsed_json
            except Exception as exc:
                parsed_chunk[key] = flat_chunk[key]
        parsed_chunks.append(parsed_chunk)

    # Finally create the martian.Records
    chunk_outs = [martian.Record(c) for c in parsed_chunks]


    # And run the stage
    env = globals()
    env.update(locals())
    exec("module.join(args, outs, chunk_defs, chunk_outs)", env, env)

    _write_outs(outs, stage.outputs)

    return outs
