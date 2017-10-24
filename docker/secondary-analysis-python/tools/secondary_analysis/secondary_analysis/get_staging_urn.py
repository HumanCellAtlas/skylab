#!/usr/bin/env python

import requests
import argparse
import time

def run(envelope_url, retry_seconds, timeout_seconds):
    start = time.time()
    urn = None
    while True:
        envelope_js = get_envelope_json(envelope_url)
        urn = get_staging_urn(envelope_js)
        if urn:
            break
        time.sleep(retry_seconds)
        current = time.time()
        if current - start >= timeout_seconds:
            message = 'Timed out while trying to get urn. Timeout seconds: {}'.format(timeout_seconds)
            raise ValueError(message)
    print(urn)

def get_envelope_json(envelope_url):
    response = requests.get(envelope_url)
    envelope_js = response.json()
    return envelope_js

def get_staging_urn(envelope_js):
    details = envelope_js.get('stagingDetails')
    if not details:
        return None
    location = details.get('stagingAreaLocation')
    if not location:
        return None
    urn = location.get('value')
    return urn

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--envelope_url', required=True)
    parser.add_argument('--retry_seconds', type=int, default=10)
    parser.add_argument('--timeout_seconds', type=int, default=600)
    args = parser.parse_args()
    run(args.envelope_url, args.retry_seconds, args.timeout_seconds)

if __name__ == '__main__':
    main()
