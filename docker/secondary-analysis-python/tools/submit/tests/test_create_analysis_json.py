import unittest
import os
import sys
import json

#pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
#sys.path.insert(0, pkg_root)

import submit.create_analysis_json as aj

class TestCreateAnalysisJson(unittest.TestCase):

    def test_create_inputs(self):
        inputs_file = self.data_file('inputs.tsv') 
        inputs = aj.create_inputs(inputs_file)
        self.assertEqual(inputs[0]['name'], 'fastq_read1')
        self.assertEqual(inputs[0]['value'], 'gs://broad-dsde-mint-dev-teststorage/path/read1.fastq.gz')
        self.assertEqual(inputs[0]['checksum'], 'd0f7d08f1980f7980f')
        self.assertEqual(inputs[1]['name'], 'fastq_read2')
        self.assertEqual(inputs[1]['value'], 'gs://broad-dsde-mint-dev-teststorage/path/read2.fastq.gz')
        self.assertEqual(inputs[1]['checksum'], 'd0f7d08f1980f7980f')
        self.assertEqual(inputs[2]['name'], 'output_prefix')
        self.assertEqual(inputs[2]['value'], 'GSM1957573')
        self.assertEqual(inputs[3]['name'], 'test_int')
        self.assertEqual(inputs[3]['value'], '123')

    def test_create_outputs(self):
        outputs_file = self.data_file('outputs.txt') 
        format_map_file = self.data_file('format_map.json') 
        outputs = aj.create_outputs(outputs_file, format_map_file)
        output_lines = []
        with open(outputs_file) as f:
            for line in f:
                output_lines.append(line.strip())
        self.assertEqual(outputs[0]['file_path'], output_lines[0])
        self.assertEqual(outputs[0]['format'], 'bam')
        self.assertEqual(outputs[0]['name'], 'Aligned.sortedByCoord.out.bam')
        self.assertEqual(outputs[1]['file_path'], output_lines[1])
        self.assertEqual(outputs[1]['format'], 'metrics')
        self.assertEqual(outputs[1]['name'], 'GSM1957573_rna_metrics')

    def test_create_outputs(self):
        outputs = aj.create_outputs(self.data_file('outputs.txt'), self.data_file('format_map.json'))

    def test_get_input_bundles(self):
        bundles = aj.get_input_bundles('foo,bar,baz')
        self.assertEqual(bundles, ['foo', 'bar', 'baz'])

    def test_get_start_end(self):
        with open(self.data_file('metadata.json')) as f:
            metadata = json.load(f)
            start, end = aj.get_start_end(metadata)
            self.assertEqual(start, '2017-09-14T19:54:11.470Z')
            self.assertEqual(end, '2017-09-14T19:54:31.871Z')

    def test_get_tasks(self):
        with open(self.data_file('metadata.json')) as f:
            metadata = json.load(f)
            tasks = aj.get_tasks(metadata)
            self.assertEqual(len(tasks), 5)
            first_task = tasks[0]
            self.assertEqual(first_task['name'], 'CollectAlignmentSummaryMetrics')
            self.assertEqual(first_task['log_out'].split('/')[-1], 'CollectAlignmentSummaryMetrics-stdout.log')
            self.assertEqual(first_task['log_err'].split('/')[-1], 'CollectAlignmentSummaryMetrics-stderr.log')
            self.assertEqual(first_task['start_time'], '2017-09-14T19:54:22.691Z')
            self.assertEqual(first_task['stop_time'], '2017-09-14T19:54:31.473Z')
            self.assertEqual(first_task['memory'], '10 GB')
            self.assertEqual(first_task['zone'], 'us-central1-b')
            self.assertEqual(first_task['cpus'], 1)
            self.assertEqual(first_task['disk_size'], 'local-disk 10 HDD')
            self.assertEqual(first_task['docker_image'], 'humancellatlas/picard')

    def test_get_format(self):
        self.assertEqual(aj.get_format('asdf', {}), 'unknown')
        self.assertEqual(aj.get_format('asdf.bam', {'.bam': 'bam'}), 'bam')
        self.assertEqual(aj.get_format('asdf.txt', {'.bam': 'bam'}), 'unknown')
        self.assertEqual(aj.get_format('asdf.bam', {'.bam': 'bam', '_metrics': 'metrics'}), 'bam')
        self.assertEqual(aj.get_format('asdf.foo_metrics', {'.bam': 'bam', '_metrics': 'metrics'}), 'metrics')

    def data_file(self, file_name):
            return os.path.split(__file__)[0] + '/data/'  + file_name

if __name__ == '__main__':
    unittest.main()
