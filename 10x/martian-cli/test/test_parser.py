"""Test the mro_parser"""
import os
import unittest

from martian_cli import mro_parser

TEST_MRO = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "test.mro")

class TestMRO(unittest.TestCase):
    """Tests against the _bcsorter_stages.mro"""

    def setUp(self):
        """Parse the MRO before each test."""
        self.mro_stages = mro_parser.get_stages(TEST_MRO)

    def test_stage_count(self):
        """All stages were parsed."""
        self.assertEqual(len(self.mro_stages), 2)

    def test_fields_count(self):
        """All IO fields are parsed."""
        empty_split_stage = self.mro_stages['EMPTY_SPLIT_STAGE']
        missing_name_stage = self.mro_stages['MISSING_NAME_STAGE']

        self.assertEqual(len(empty_split_stage.inputs), 3)
        self.assertEqual(len(empty_split_stage.outputs), 2)
        self.assertEqual(len(empty_split_stage.splits), 0)

        self.assertEqual(len(missing_name_stage.inputs), 1)
        self.assertEqual(len(missing_name_stage.outputs), 1)
        self.assertEqual(len(missing_name_stage.splits), 1)

    def test_default_field(self):
        """Nameless field parsed correctly."""
        missing_name_stage = self.mro_stages['MISSING_NAME_STAGE']

        self.assertEqual(missing_name_stage.inputs[0].name, "default")
        self.assertEqual(missing_name_stage.outputs[0].name, "default")
        self.assertEqual(missing_name_stage.splits[0].name, "default")

    def test_types(self):
        """Types are parsed correctly."""
        empty_split_stage = self.mro_stages['EMPTY_SPLIT_STAGE']
        int_input = [i for i in empty_split_stage.inputs
                     if i.name == "int_input"][0]
        fastq_input = [i for i in empty_split_stage.inputs
                       if i.name == "fastq_array_input"][0]
        int_input2 = [i for i in empty_split_stage.inputs
                      if i.name == "int_input2"][0]
        map_output = [i for i in empty_split_stage.outputs
                      if i.name == "map_array_output"][0]
        float_output = [i for i in empty_split_stage.outputs
                        if i.name == "float_output"][0]

        self.assertEqual(int_input.type, "int")
        self.assertEqual(fastq_input.type, "fastq[]")
        self.assertEqual(int_input2.type, "int")
        self.assertEqual(map_output.type, "map[]")
        self.assertEqual(float_output.type, "float")
