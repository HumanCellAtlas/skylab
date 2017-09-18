#!/usr/bin/env python

import requests
import json
import argparse

def run(submit_url, analysis_json_path):
    print('Starting submission with url: {0}'.format(submit_url))

    # 1. Get envelope url
    response = requests.get(submit_url)
    envelope_url = get_entity_url(response.json(), 'submissionEnvelopes')

    # 2. Create envelope, get analysis and submission urls
    print('Creating submission envelope at {0}'.format(envelope_url))
    response = requests.post(envelope_url, '{}')
    envelope_js = response.json()
    analyses_url = get_entity_url(envelope_js, 'analyses')
    submission_url = get_entity_url(envelope_js, 'submissionEnvelope')

    # 3. Create analysis, get input bundles url, file refs url
    print('Creating analysis at {0}'.format(analyses_url))
    json_header = {'Content-type': 'application/json'}
    with open(analysis_json_path) as f:
        analysis_json_contents = json.load(f)
    response = requests.post(analyses_url, headers = json_header, data = analysis_json_contents)
    print(response.text)
    analysis_js = response.json()
    input_bundles_url = get_entity_url(analysis_js, 'add-input-bundles')
    file_refs_url = get_entity_url(analysis_js, 'add-file-reference')

    # 4. Add input bundles
    print('Adding input bundles at {0}'.format(input_bundles_url))
    input_bundle_uuid = get_input_bundle_uuid(analysis_json_contents)
    response = requests.post(input_bundles_url, data = {"bundleUuids": [input_bundle_uuid]})
    print(response.text)

    # 5. Add file references
    print('Adding file references at {0}'.format(file_refs_url))
    output_files = get_output_files(analysis_json_contents)
    response = requests.put(file_refs_url, headers = json_header, data = output_files)
    print(response.text)

    # 6. Confirm submission
    confirmation_url = '{0}/{1}'.format(submission_url, 'confirmation')
    print('Confirming submission at {0}'.format(confirmation_url))
    response = requests.put(confirmation_url, headers = json_header)
    print(response)

def get_entity_url(js, entity):
    entity_url = js['_links'][entity]['href'].split('{')[0]
    print('Got url for {0}: {1}'.format(entity, entity_url))
    return entity_url

def get_input_bundle_uuid(analysis_json):
    bundle = analysis_json['input_bundles'][0]
    uuid = bundle.split('/')[-2]
    print('Input bundle uuid {0}'.format(uuid))
    return uuid

def get_output_files(analysis_json):
    outputs = analysis_json['outputs']
    output_refs = []
    for o in outputs:
        output_ref = {}
        file_name = o['file_path'].split('/')[-1] 
        output_ref['fileName'] = file_name
        content = {}
        content['name'] = file_name
        content['format'] = o['format']
        output_ref['content'] = content
        output_refs.append(output_ref)
    return output_refs
        
#if __name__ == '__main__':
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-submit_url')
    parser.add_argument('-analysis_json_path')
    args = parser.parse_args()
    #print(args.submit_url, args.analysis_json_path)
    run(args.submit_url, args.analysis_json_path)
