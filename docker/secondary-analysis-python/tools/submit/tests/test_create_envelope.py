import unittest
import os
import json

import submit.create_envelope as submit

class TestCreateEnvelope(unittest.TestCase):

    def test_get_entity(self):
        with open(self.data_file('response.json')) as f:
            js = json.load(f)
            entity_url = submit.get_entity_url(js, 'analyses')
            self.assertEqual(entity_url, 'http://api.ingest.dev.data.humancellatlas.org/analyses')

    def test_get_input_bundle_uuid(self):
        with open(self.data_file('analysis.json')) as f:
            js = json.load(f)
            self.assertEqual(submit.get_input_bundle_uuid(js), '75a7f618-9adc-48af-a249-0010305160f6')

    def test_get_output_files(self):
        with open(self.data_file('analysis.json')) as f:
            js = json.load(f)
            outputs = submit.get_output_files(js)
            self.assertEqual(len(outputs), 3)
            self.assertEqual(outputs[0]['fileName'], 'sample.bam')
            self.assertEqual(outputs[0]['content']['name'], 'sample.bam')
            self.assertEqual(outputs[0]['content']['format'], 'bam')

    def test_check_status_bad_codes(self):
        with self.assertRaises(ValueError):
            submit.check_status(404, 'foo')
        with self.assertRaises(ValueError):
            submit.check_status(500, 'foo')
        with self.assertRaises(ValueError):
            submit.check_status(301, 'foo')

    def test_check_status_acceptable_codes(self):
        try:
          submit.check_status(200, 'foo')
        except ValueError as e:
            self.fail(str(e))

        try:
            submit.check_status(202, 'foo')
        except ValueError as e:
            self.fail(str(e))

        try:
            submit.check_status(301, 'foo', '3xx')
        except ValueError as e:
            self.fail(str(e))

    def data_file(self, file_name):
        return os.path.split(__file__)[0] + '/data/'  + file_name

if __name__ == '__main__':
    unittest.main()
