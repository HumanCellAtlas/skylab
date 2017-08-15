"""Parses an MRO file into a python representation.

There are two pieces to this. First, pyparsing creates a grammar
for the mro format and uses that to tokenize it into
python objects. Second, those objects are tidied up into something
that's a somewhat intuitive interface.
"""

import collections
import os
import re

import pyparsing

from martian_cli import common

__all__ = ["get_stages"]

# A few little types for various Martian entities
MartianIOField = collections.namedtuple(
    "MartianIOField",
    ["modifier",     # A string in ["in", "out"]
     "type",         # A string like "json", "int", "float", "map[]", etc.
     "name",         # A string that identifies the field
     "help"])        # An optional help message string

MartianStage = collections.namedtuple(
    "MartianStage",
    ["name",         # A string that identifies the stage
     "inputs",       # A list of MartianIOFields
     "outputs",      # A list of MartianIOFields
     "splits",       # A list of MartianIOFields
     "source",       # A string that is a path to a directory that contains an __init__.py
     "force_split"]) # Run split and join even if the splits are empty

MartianCall = collections.namedtuple(
    "MartianCall",
    ["name",         # A string identifying callable that is being called
     "assignments"]) # A dict with keys for input fields of the callable and values for the
                     #   assigned value. Note that the value needs some parsing, so "self.lanes"
                     #   refers to the "lanes" input to the pipeline, etc.

MartianPipeline = collections.namedtuple(
    "MartianPipeline",
    ["name",         # A string identifying the pipeline
     "inputs",       # A list of MartianIOFields
     "outputs",      # A list of MartianIOFields
     "calls",        # A list of MartianCalls
     "return_"])     # A dict with keys for the pipeline's output names and values for the
                     #   assignment of call outputs to pipeline outputs


def get_stages(mro_file_path):
    """Get MartianStages from an MRO file.

    Inputs:
        mro_file_path: Path to an MRO file that defines a pipeline

    Returns:
        A dict where keys are stage names and values are MartianStages
    """

    stages = {}
    parsed_mro = _mro_parser().parseFile(mro_file_path)

    if parsed_mro.stages:
        stages_ = _create_stages(parsed_mro, mro_file_path)
        for stage in stages_:
            stages[stage.name] = stage

    return stages


def _create_stages(parsed_mro, mro_file_path):
    """Create a list of MartianStages from parsed MRO.

    Args:
      parsed_mro - The outputs of a parse_mro call
      mro_file_path - The path to the MRO file that was parsed. This is needed to get a full path
        to the stage's source.

    Returns:
      A list of MartianStages from the parsed MRO file.
    """

    stages = []

    for parsed_stage in parsed_mro.stages:
        stage_name = parsed_stage[0]

        def normalize_entry(entry):
            """IO entries in Martian stages can have comments or not. *Intriguingly*, they can
            also be anonymous. This handles that variation and stick in empty values things
            that are missing.
            """
            if len(entry) == 3:
                return list(entry) + [None]
            elif len(entry) == 2:  # Looking at you, SORT_BY_BC
                return list(entry) + ["default", None]
            else:
                return list(entry)

        normalized_entries = [normalize_entry(e) for e in parsed_stage.stage_entries]
        normalized_splits = [normalize_entry(e) for e in parsed_stage.split]

        stage_inputs = [MartianIOField(*i) for i in normalized_entries if i[0] == 'in']
        stage_outputs = [MartianIOField(*i) for i in normalized_entries if i[0] == 'out']
        stage_splits = [MartianIOField(*i) for i in normalized_splits]
        stage_source = [os.path.join(os.path.dirname(mro_file_path), re.sub(r'^"|"$', "", i[2]))
                        for i in parsed_stage.stage_entries if i[0] == 'src'][0]

        stages.append(
            MartianStage(
                stage_name,
                stage_inputs,
                stage_outputs,
                stage_splits,
                stage_source,
                parsed_stage.split != ''))

    return stages


def _mro_parser():
    """Create and return a pyparsing parser for mro files."""

    # A few helpful pyparsing constants
    EQUALS, SEMI, LBRACE, RBRACE, LPAREN, RPAREN = map(pyparsing.Suppress, '=;{}()')
    mro_label = pyparsing.Word(pyparsing.alphanums + '_')
    mro_modifier = pyparsing.oneOf(["in", "out", "src"])
    mro_type = pyparsing.oneOf([
        "bool", "bool[]",
        "int", "int[]",
        "float", "float[]",
        "map", "map[]",
        "string", "string[]", "string[][]",
        "path", "path[]",
        "py"] +
        common.MARTIAN_FILETYPES + [x + '[]' for x in common.MARTIAN_FILETYPES])

    # First parse includes
    include = pyparsing.Literal("@include").suppress() + pyparsing.quotedString
    includes = pyparsing.ZeroOrMore(include).setResultsName("includes")
    includes.addParseAction(pyparsing.removeQuotes)

    # Then parse filetypes
    filetype = (pyparsing.Literal("filetype").suppress() +
                pyparsing.oneOf(common.MARTIAN_FILETYPES) +
                SEMI)
    filetypes = pyparsing.ZeroOrMore(filetype).setResultsName("filetypes")

    #
    # Stage
    #

    # Now define the parts of a stage

    # First we have a "stage entry", which is a line in the stage body, it looks like "in int lane"
    stage_entry = pyparsing.Group(
        mro_modifier +
        mro_type +
        pyparsing.Optional(pyparsing.Word(pyparsing.printables, excludeChars=',')) +
        pyparsing.Optional(pyparsing.QuotedString('"')))

    # Note that stage entries a comma-delimited, but there's a trailing comma so we need the
    # pyparsing.Empty option for matching
    stage_entries = pyparsing.delimitedList(pyparsing.Or([stage_entry, pyparsing.Empty()]))

    # Each stage can have two parts, the main part and a "split using" part
    split = (pyparsing.Literal("split using").suppress() +
             LPAREN +
             pyparsing.Optional(pyparsing.Group(stage_entries).setResultsName("split")) +
             RPAREN)

    stage = pyparsing.Group(
        pyparsing.Literal("stage").suppress() +
        mro_label +
        LPAREN +
        pyparsing.Group(stage_entries).setResultsName("stage_entries") +
        RPAREN +
        pyparsing.Optional(split))

    # Now create a dict of the stages, with the MRO labels for keys
    stages = pyparsing.Dict(pyparsing.ZeroOrMore(stage)).setResultsName("stages")

    #
    # Pipeline
    #

    # Calls
    call_entry = pyparsing.Group(
        pyparsing.Word(pyparsing.printables, excludeChars="=") +
        EQUALS +
        pyparsing.Word(pyparsing.printables, excludeChars=','))

    call_entries = pyparsing.delimitedList(pyparsing.Or([call_entry, pyparsing.Empty()]))

    call_modifier = pyparsing.oneOf(["local", "preflight"])

    call = pyparsing.Group(
        pyparsing.Literal("call").suppress() +
        pyparsing.ZeroOrMore(call_modifier).suppress() +
        mro_label +
        LPAREN +
        pyparsing.Group(call_entries).setResultsName("call_entries") +
        RPAREN)

    calls = pyparsing.Dict(pyparsing.ZeroOrMore(call)).setResultsName("pipeline_calls")

    # Return
    return_entry = call_entry
    return_entries = pyparsing.delimitedList(pyparsing.Or([return_entry, pyparsing.Empty()]))
    return_ = (pyparsing.Literal("return").suppress() +
               LPAREN +
               pyparsing.Group(return_entries).setResultsName("pipeline_return") +
               RPAREN)

    # Pipeline header
    pipeline_header_entry = pyparsing.Group(
        mro_modifier +
        mro_type +
        pyparsing.Word(pyparsing.printables, excludeChars=",") +
        pyparsing.Optional(pyparsing.quotedString))

    pipeline_header_entries = pyparsing.delimitedList(
        pyparsing.Or([pipeline_header_entry, pyparsing.Empty()]))

    pipeline = (pyparsing.Literal("pipeline").suppress() +
                mro_label.setResultsName("pipeline_name") +
                LPAREN +
                pyparsing.Group(pipeline_header_entries).setResultsName("pipeline_header") +
                RPAREN +
                LBRACE +
                calls +
                return_ +
                RBRACE)

    mro_file = pyparsing.Each(
        [pyparsing.Optional(includes), filetypes, stages, pyparsing.Optional(pipeline)])
    mro_file.ignore(pyparsing.pythonStyleComment)

    return mro_file
