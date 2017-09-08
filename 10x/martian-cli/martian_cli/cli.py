"""Module to define and expose the CLI to Martian

Contains the entry point for the package. Provides subcommands for running
martian stages as well as looking at what stages exists and what their inputs
and ouptuts are.
"""

import argparse
import json
import os
import re
import sys
import traceback

from martian_cli import adapters, common, mro_parser


def verify_environment():
    """Check that required environment values are present.

    This is just checking for MROPATH and the ability to import
    martian.
    """

    if "MROPATH" not in os.environ:
        raise EnvironmentError(
            "MROPATH is not in the environment. You probably need to source "
            "sourceme.bash in cellranger before running this tool.")

    try:
        import martian
    except ImportError:
        print sys.path
        traceback.print_exc()
        raise ImportError(
            "Could not import martian. You probably need to source "
            "sourceme.bash in cellranger before running this tool.")


def stage_list(args):
    """Function that handles `martian stage list`.

    Just prints out the loaded stages.
    """

    for stage in args.stages:
        print stage


def stage_describe(args):
    """Function that handles `martian stage describe <stage_name>`

    Prints out inputs, outputs, and splits for the stage.
    """

    stage = args.stages[args.stage_name]

    print "Stage:", stage.name, "\n"
    print "Inputs:"
    for input_ in stage.inputs:
        print "\t", input_.name, input_.type, input_.help or ""

    print "\n"
    print "Outputs:"
    for output in stage.outputs:
        print "\t", output.name, output.type, output.help or ""

    print "\n"
    print "Splits:"
    for split in stage.splits:
        print "\t", split.name, split.type, split.help or ""

    print "\n"
    print "Source directory:"
    print "\t", stage.source


def stage_run(args):
    """Function that handles `martian stage run <stage> <phase> ...`

    Generally just passes everything to the relevant function in
    adapters.
    """

    print "stage_run args:", args
    if args.phase == "split":
        adapters.split(args.stage, args)
    if args.phase == "main":
        adapters.main(args.stage, args)
    if args.phase == "join":

        #for key in args.__dict__:
        #    if isinstance(getattr(args, key), list) and len(getattr(args, key)) == 1:
        #        setattr(args, key, getattr(args, key)[0])

        adapters.join(args.stage, args)


def get_parser(stages):
    """Create the argument parser for the CLI.

    Args:
      stages - a list of all MartianStages in MRO files in the MROPATH
    """

    # martian
    parser = argparse.ArgumentParser(prog="martian")
    subparsers = parser.add_subparsers()

    # martian stage
    stage_parser = subparsers.add_parser(
        "stage", help="Work with Martian stages.")

    stage_subparsers = stage_parser.add_subparsers(
        title="Stage subcommands",
        help="Actions than can be performed on Martian stages.")

    # martian stage list
    stage_list_parser = stage_subparsers.add_parser(
        "list",
        help="List all available stages.")
    stage_list_parser.set_defaults(func=stage_list)
    stage_list_parser.set_defaults(stages=stages)

    # martian stage describe <stage_name>
    stage_describe_parser = stage_subparsers.add_parser(
        "describe",
        help="Describe the inputs, outputs, and source location of a stage")
    stage_describe_parser.add_argument(
        "stage_name",
        help="Name of the stage to describe")
    stage_describe_parser.set_defaults(func=stage_describe)
    stage_describe_parser.set_defaults(stages=stages)

    # martian stage run <stage_name> <stage_phase> <stage_args...>
    stage_run_parser = stage_subparsers.add_parser(
        "run",
        help="Run a stage")

    stage_run_subparsers = stage_run_parser.add_subparsers(
        title="Stages available to be run",
        help="Names of available stages.")

    for stage in stages.values():

        individual_stage_parser = stage_run_subparsers.add_parser(
            stage.name,
            help="Execute stage " + stage.name)
        individual_stage_subparsers = individual_stage_parser.add_subparsers()

        # Some stages don't have a split or join
        available_stage_phases = ['split', 'join', 'main'] if (stage.splits or stage.force_split) else ['main']

        for phase in available_stage_phases:
            phase_parser = individual_stage_subparsers.add_parser(
                phase,
                help='Run the ' + phase + ' of ' + stage.name)

            phase_parser.set_defaults(func=stage_run)
            phase_parser.set_defaults(phase=phase)
            phase_parser.set_defaults(stage=stage)

            for input_ in _stage_inputs(stage, phase):
                help_message = "Type: " + input_.type
                if input_.help:
                    help_message += " Help: " + input_.help
                phase_parser.add_argument(
                    "--" + input_.name,
                    type=martian_type_to_python_type(input_.type),
                    nargs=martian_type_to_nargs(input_.type),
                    default=None,
                    help=help_message)

            # Handle the "split_file" for mains that come after a split
            if phase == 'main' and 'split' in available_stage_phases:
                phase_parser.add_argument(
                    "--split_file",
                    type=martian_type_to_python_type("File"),
                    nargs=martian_type_to_nargs("File"),
                    default=None,
                    help="File with split arguments.")

    return parser


def _strip_quotes(file_arg):
    """Strip leading and trailing quotes.

    These are showing up in paths...
    """
    return re.sub("^[\'\"]|[\'\"]$", "", file_arg)

def _load_map(arg):
    """Try to load a map argument.

    This seems to show up in a couple of different ways in WDL. It can just be a json string, or
    sometimes it get's written to a 2 column tsv.
    """

    if os.path.isfile(arg):
        output = {}
        with open(arg) as arg_f:
            for line in arg_f:
                key, val = line.strip().split('\t')
                output[key] = val
        return output
    return json.loads(arg)

def _str_to_bool(val):
    """Convert a bool-like string to a bool. Useful for exposing options in the cli."""
    if val.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif val.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def _stage_inputs(stage, phase):
    """Get the inputs to a phase of a stage.

    Each phase (split, main, join) of a stage takes different input fields, and the types
    of the fields can also change based on phase.
    """

    def arrayify(martian_io_field):
        """Convert the type of a Martian input field to an array of that type.

        This is necessary for the join phase.
        """
        return mro_parser.MartianIOField(
            martian_io_field.modifier,
            martian_io_field.type + '[]',
            martian_io_field.name,
            martian_io_field.help)

    def add_tag_to_name(martian_io_field, tag):
        return mro_parser.MartianIOField(
            martian_io_field.modifier,
            martian_io_field.type,
            martian_io_field.name + '_' + tag,
            martian_io_field.help)

    if phase == 'split':
        return stage.inputs
    elif phase == 'main':
        return stage.inputs + stage.splits
    elif phase == 'join':
        # The inputs to join are arrays of the split and output fields since it's pulling
        # together outputs of multiple main steps.
        # Also, "split" and "output" need to be added to the field names or there are collisions
        return stage.inputs + \
            [add_tag_to_name(arrayify(s), "split") for s in stage.splits] + \
            [add_tag_to_name(arrayify(s), "output") for s in stage.outputs]


def martian_type_to_python_type(martian_type):
    """Convert a Martian type to a python type for argparse."""
    martian_to_python = {"path": _strip_quotes,
                         "int": int,
                         "bool": _str_to_bool,
                         "float": float,
                         "string": str,
                         "File": _strip_quotes,
                         "map": _load_map}
    for ftype in common.MARTIAN_FILETYPES:
        martian_to_python[ftype] = str

    return martian_to_python[martian_type.strip('[]')]


def martian_type_to_nargs(martian_type):
    """Get the argparse nargs parameter for a Martian type."""

    if martian_type.endswith("[]"):
        return "*"
    else:
        return '?'


def load_stages():
    """Load all stages defined in MRO files in the MROPATH."""

    def load_stages_from_dir(mro_dir):
        """Iterate over MRO file in a directory, parse them, and accumulate
        their stages.
        """
        stages = {}
        for file_name in os.listdir(mro_dir):
            if file_name.endswith(".mro"):
                stages.update(mro_parser.get_stages(
                    os.path.join(mro_dir, file_name)))
        return stages

    stages = {}
    for mro_dir in os.environ["MROPATH"].split(':'):
        stages.update(load_stages_from_dir(mro_dir))
    return stages


def main():
    """Entry point for the package."""

    # Make sure we're running in the cellranger environment
    verify_environment()

    # Load all the available stages
    stages = load_stages()

    # Create an argument parser out of the stages
    parser = get_parser(stages)

    args = parser.parse_args()
    args.func(args)
