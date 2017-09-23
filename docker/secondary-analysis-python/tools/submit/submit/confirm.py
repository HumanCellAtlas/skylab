#!/usr/bin/env python

import requests
import json
import argparse
import stage
import time

def run(envelope_url):
    start = time.time()
    current = start
    while current - start < 120:
        print('Getting status for {}'.format(envelope_url))
        envelope_js = get_envelope_json(envelope_url)
        status = envelope_js.get('submissionState')
        print('submissionState: {}'.format(status))
        if status == 'Valid':
            break
        time.sleep(10)
        current = time.time()
    print('Confirming submission')
    headers = {
        'Content-type': 'application/json'
    }
    response = requests.put('{}/submissionEvent'.format(envelope_url), headers = headers)
    print(response.text)

def get_envelope_json(envelope_url):
    response = requests.get(envelope_url)
    envelope_js = response.json()
    return envelope_js

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-envelope_url')
    args = parser.parse_args()
    run(args.envelope_url)

if __name__ == '__main__':
    main()
