import unittest
import os
import sys
import json
import submit.get_staging_urn as gsu

class TestGetStagingUrn(unittest.TestCase):

    def test_empty_js(self):
        js = {}
        self.assertIsNone(gsu.get_staging_urn(js))

    def test_null_details(self):
        js = { 
            'stagingDetails': None
        }
        self.assertIsNone(gsu.get_staging_urn(js))

    def test_null_location(self):
        js = { 
            'stagingDetails': {
                'stagingAreaLocation': None
            }
        }
        self.assertIsNone(gsu.get_staging_urn(js))

    def test_null_value(self):
        js = { 
            'stagingDetails': {
                'stagingAreaLocation': {
                    'value': None
                }
            }
        }
        self.assertIsNone(gsu.get_staging_urn(js))

    def test_valid_value(self):
        js = { 
            'stagingDetails': {
                'stagingAreaLocation': {
                    'value': 'test_urn'
                }
            }
        }
        self.assertEqual(gsu.get_staging_urn(js), 'test_urn')

if __name__ == '__main__':
    unittest.main()
