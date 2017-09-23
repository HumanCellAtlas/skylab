#!/usr/bin/env python

import requests
import json
import argparse
import stage
import time

def run(envelope_url):
    start = time.time()
    current = start
    urn = None
    while current - start < 120:
        envelope_js = get_envelope_json(envelope_url)
        urn = get_submission_urn(envelope_js)
        if urn:
            break
        time.sleep(10)
        current = time.time()
    print(urn)

def get_envelope_json(envelope_url):
    response = requests.get(envelope_url)
    envelope_js = response.json()
    return envelope_js

def get_submission_urn(envelope_js):
    details = envelope_js.get('stagingDetails')
    if not details:
        return None
    location = details.get('stagingAreaLocation')
    if not location:
        return None
    urn = location.get('value')
    if not urn:
        return None
    return urn

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-envelope_url')
    args = parser.parse_args()
    run(args.envelope_url)

if __name__ == '__main__':
    main()
