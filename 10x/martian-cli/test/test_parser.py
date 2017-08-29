"""Test the mro_parser"""
import os
import unittest

from martian_cli import mro_parser

BCSORTER_MRO = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "_bcsorter_stages.mro")
MAKE_FASTQS_MRO = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "_make_fastqs_stages.mro")

class TestBCSorterMRO(unittest.TestCase):
    """Tests against the _bcsorter_stages.mro"""

    def setUp(self):
        """Parse the MRO before each test."""
        self.mro_stages = mro_parser.get_stages(BCSORTER_MRO)

    def test_stage_count(self):
        """All stages were parsed."""
        self.assertEqual(len(self.mro_stages), 2)

    def test_fields_count(self):
        """All IO fields are parsed."""
        sort_by_bc_stage = self.mro_stages['SORT_BY_BC']
        bucket_by_bc_stage = self.mro_stages['BUCKET_BY_BC']

        self.assertEqual(len(sort_by_bc_stage.inputs), 1)
        self.assertEqual(len(sort_by_bc_stage.outputs), 2)
        self.assertEqual(len(sort_by_bc_stage.splits), 2)

        self.assertEqual(len(bucket_by_bc_stage.inputs), 3)
        self.assertEqual(len(bucket_by_bc_stage.outputs), 1)
        self.assertEqual(len(bucket_by_bc_stage.splits), 2)

    def test_default_field(self):
        """Nameless field parsed correctly."""
        sort_by_bc_stage = self.mro_stages['SORT_BY_BC']

        bam_output = [i for i in sort_by_bc_stage.outputs if i.type == "bam"][0]

        self.assertEqual(bam_output.name, "default")

class TestMakeFastqsMRO(unittest.TestCase):
    """Tests against the _make_fastqs_stage.mro"""

    def setUp(self):
        """Parse the MRO before each test."""
        self.mro_stages = mro_parser.get_stages(MAKE_FASTQS_MRO)

    def test_stage_count(self):
        """All stages were parsed."""
        self.assertEqual(len(self.mro_stages), 6)

    def test_fields_count(self):
        """All IO fields are parsed."""
        make_fastqs_preflight_local_stage = self.mro_stages["MAKE_FASTQS_PREFLIGHT_LOCAL"]
        make_fastqs_preflight_stage = self.mro_stages["MAKE_FASTQS_PREFLIGHT"]
        prepare_samplesheet = self.mro_stages["PREPARE_SAMPLESHEET"]


        self.assertEqual(len(make_fastqs_preflight_local_stage.inputs), 15)
        self.assertEqual(len(make_fastqs_preflight_local_stage.outputs), 0)
        self.assertEqual(len(make_fastqs_preflight_local_stage.splits), 0)

        self.assertEqual(len(make_fastqs_preflight_stage.inputs), 6)
        self.assertEqual(len(make_fastqs_preflight_stage.outputs), 0)
        self.assertEqual(len(make_fastqs_preflight_stage.splits), 0)

        self.assertEqual(len(prepare_samplesheet.inputs), 6)
        self.assertEqual(len(prepare_samplesheet.outputs), 3)
        self.assertEqual(len(prepare_samplesheet.splits), 0)

    def test_types(self):
        """Types are parsed correctly."""
        merge_fastqs_by_lane_sample_stage = self.mro_stages["MERGE_FASTQS_BY_LANE_SAMPLE"]

        fastq_path_input = [i for i in merge_fastqs_by_lane_sample_stage.inputs
                            if i.name == "fastq_path"][0]
        files_merged_output = [i for i in merge_fastqs_by_lane_sample_stage.outputs
                               if i.name == "files_merged"][0]
        merged_file_paths_output = [i for i in merge_fastqs_by_lane_sample_stage.outputs
                                    if i.name == "merged_file_paths"][0]
        lane_split = [i for i in merge_fastqs_by_lane_sample_stage.splits
                      if i.name == "lane"][0]

        self.assertEqual(fastq_path_input.type, "path")
        self.assertEqual(files_merged_output.type, "bool")
        self.assertEqual(merged_file_paths_output.type, "string[]")
        self.assertEqual(lane_split.type, "int")
