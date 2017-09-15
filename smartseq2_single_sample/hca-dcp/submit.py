#!/usr/bin/env python

import requests
import json

def run(submit_url):
    # Get envelope url
    response = requests.get(submit_url)
    js = response.json()
    envelope_url = js['_links']['submissionEnvelopes']['href'].split('{')[0]
    print(envelope_url)

    # Create envelope, get analysis and submission urls
    response = requests.post(envelope_url, '{}')
    envelope_js = response.json()
    analyses_url = envelope_js['_links']['analyses']['href'].split('{')[0]
    submission_url = envelope_js['_links']['submissionEnvelope']['href'].split('{')[0]
    print(response.text)
    print(submission_url)
    print(analyses_url)

    # Create analysis, get input bundles url, file refs url
    headers = {
      'Content-type': 'application/json'
    }
    response = requests.post(analyses_url, headers = headers, data = '{}')
    analysis_js = response.json()
    input_bundles_url = analysis_js['_links']['add-input-bundles']['href'].split('{')[0]
    print(input_bundles_url)
    file_refs_url = analysis_js['_links']['add-file-reference']['href'].split('{')[0]
    print(file_refs_url)

    # Add input bundles
    response = requests.post(input_bundles_url, data = {"bundleUuids": ["46cf0282-ac03-4c57-9274-c8e1bdd2aed7"]})
    print(response.text)

    # Add file references
    file_ref_d = {
        "fileName": "ERR1630013.fastq.gz",
        "content": {
            "lane": 1,
            "type": "reads",
            "name": "ERR1630013.fastq.gz",
            "format": ".fastq.gz"
        }
    }
    #response = requests.put(file_refs_url, headers = headers, data = file_ref_d)
    #print(response.text)

    # Confirm submission
    confirmation_url = '{0}/{1}'.format(submission_url, 'confirmation')
    print(confirmation_url)
    response = requests.put(confirmation_url, headers = headers)
    print(response)

if __name__ == '__main__':
    run('http://api.ingest.dev.data.humancellatlas.org/')
