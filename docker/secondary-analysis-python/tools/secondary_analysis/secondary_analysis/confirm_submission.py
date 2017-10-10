#!/usr/bin/env python

import requests
import argparse
import time

def run(envelope_url, retry_seconds, timeout_seconds):
    start = time.time()
    while True:
        print('Getting status for {}'.format(envelope_url))
        envelope_js = get_envelope_json(envelope_url)
        status = envelope_js.get('submissionState')
        print('submissionState: {}'.format(status))
        if status == 'Valid':
            break
        time.sleep(retry_seconds)
        current = time.time()
        if current - start >= timeout_seconds:
            message = 'Timed out while waiting for Valid status. Timeout seconds: {}'.format(timeout_seconds)
            raise ValueError(message)
    confirm(envelope_url)

def confirm(envelope_url):
    print('Confirming submission')
    headers = {
        'Content-type': 'application/json'
    }
    response = requests.put('{}/submissionEvent'.format(envelope_url), headers=headers)
    print(response.text)

def get_envelope_json(envelope_url):
    response = requests.get(envelope_url)
    envelope_js = response.json()
    return envelope_js

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--envelope_url', required=True)
    parser.add_argument('--retry_seconds', type=int, default=10)
    parser.add_argument('--timeout_seconds', type=int, default=600)
    args = parser.parse_args()
    run(args.envelope_url, args.retry_seconds, args.timeout_seconds)

if __name__ == '__main__':
    main()
