import unittest
import os
import sys
import json

#pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
#sys.path.insert(0, pkg_root)

import submit.submit as submit

class TestService(unittest.TestCase):

    def test_get_entity(self):
        with open(self.data_file('response.json')) as f:
            js = json.load(f)
            entity_url = submit.get_entity_url(js, 'analyses')
            self.assertEqual(entity_url, 'http://api.ingest.dev.data.humancellatlas.org/analyses')

    def test_get_input_bundle_uuid(self):
        with open(self.data_file('analysis.json')) as f:
            js = json.load(f)
            self.assertEqual(submit.get_input_bundle_uuid(js), '23bd7eb0-23a3-4898-b043-f7e982de281f')

    def test_get_output_files(self):
        with open(self.data_file('analysis.json')) as f:
            js = json.load(f)
            outputs = submit.get_output_files(js)
            self.assertEqual(len(outputs), 3)
            self.assertEqual(outputs[0]['fileName'], 'sample.bam')
            self.assertEqual(outputs[0]['content']['name'], 'sample.bam')
            self.assertEqual(outputs[0]['content']['format'], 'bam')

    def data_file(self, file_name):
        return os.path.split(__file__)[0] + '/data/'  + file_name

if __name__ == '__main__':
    unittest.main()
