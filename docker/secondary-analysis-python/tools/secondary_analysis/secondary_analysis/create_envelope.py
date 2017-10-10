#!/usr/bin/env python

import requests
import json
import argparse

def run(submit_url, analysis_json_path):
    # 1. Get envelope url
    print('Getting envelope url from {}'.format(submit_url))
    response = requests.get(submit_url)
    check_status(response.status_code, response.text)
    envelope_url = get_entity_url(response.json(), 'submissionEnvelopes')

    # 2. Create envelope, get analysis and submission urls
    print('Creating submission envelope at {0}'.format(envelope_url))
    response = requests.post(envelope_url, '{}')
    check_status(response.status_code, response.text)
    envelope_js = response.json()
    analyses_url = get_entity_url(envelope_js, 'analyses')
    print('Creating analysis at {0}'.format(analyses_url))
    submission_url = get_entity_url(envelope_js, 'submissionEnvelope')
    with open('submission_url.txt', 'w') as f:
        f.write(submission_url)

    # 3. Create analysis, get input bundles url, file refs url
    json_header = {'Content-type': 'application/json'}
    with open(analysis_json_path) as f:
        analysis_json_contents = json.load(f)
    response = requests.post(analyses_url, headers=json_header, data=json.dumps(analysis_json_contents))
    check_status(response.status_code, response.text)
    analysis_js = response.json()
    input_bundles_url = get_entity_url(analysis_js, 'add-input-bundles')
    file_refs_url = get_entity_url(analysis_js, 'add-file-reference')

    # 4. Add input bundles
    print('Adding input bundles at {0}'.format(input_bundles_url))
    input_bundle_uuid = get_input_bundle_uuid(analysis_json_contents)
    bundle_refs_js = json.dumps({"bundleUuids": [input_bundle_uuid]}, indent=2)
    print(bundle_refs_js)
    response = requests.put(input_bundles_url, headers=json_header, data=bundle_refs_js)
    check_status(response.status_code, response.text)

    # 5. Add file references
    print('Adding file references at {0}'.format(file_refs_url))
    output_files = get_output_files(analysis_json_contents)
    for file_ref in output_files:
        print('Adding file: {}'.format(file_ref['fileName']))
        response = requests.put(file_refs_url, headers=json_header, data=json.dumps(file_ref))
        check_status(response.status_code, response.text)

def check_status(status, response_text, expected='2xx'):
    """Check that status is in range 200-299 or the specified range, if given.
    Raises a ValueError and prints response_text if status is not in the expected range. Otherwise,
    just returns silently.
    Args:
        status (int): The actual HTTP status code.
        response_text (str): Text to print along with status code when mismatch occurs
        expected (str): The range of acceptable values represented as a string.
    Examples:
        check_status(200, 'foo') passes
        check_status(404, 'foo') raises error
        check_status(301, 'bar') raises error
        check_status(301, 'bar', '3xx') passes
    """
    first_digit = int(expected[0])
    low = first_digit * 100
    high = low + 99
    matches = low <= status <= high
    if not matches:
        message = 'HTTP status code {0} is not in expected range {1}. Response: {2}'.format(status, expected, response_text)
        raise ValueError(message)

def get_entity_url(js, entity):
    entity_url = js['_links'][entity]['href'].split('{')[0]
    print('Got url for {0}: {1}'.format(entity, entity_url))
    return entity_url

def get_input_bundle_uuid(analysis_json):
    bundle = analysis_json['input_bundles'][0]
    uuid = bundle
    print('Input bundle uuid {0}'.format(uuid))
    return uuid

def get_output_files(analysis_json):
    outputs = analysis_json['outputs']
    output_refs = []
    for o in outputs:
        output_ref = {}
        file_name = o['file_path'].split('/')[-1] 
        output_ref['fileName'] = file_name
        output_ref['content'] = {
            'name': file_name,
            'format': o['format']
        }
        output_refs.append(output_ref)
    return output_refs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--submit_url', required=True)
    parser.add_argument('--analysis_json_path', required=True)
    args = parser.parse_args()
    run(args.submit_url, args.analysis_json_path)

if __name__ == '__main__':
    main()
