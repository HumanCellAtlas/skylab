"""Runs a single portability test from a test directory."""

import argparse
import datetime
import glob
import json
import os
import sys
import time

import requests

# This is what we need to find in a test directory. I guess we'll glob the PR
# test WDL itself since it's named after the pipeline.
TEST_DIR_LAYOUT = {
    "inputs": "test_inputs.json",
    "dependencies": "dependencies.json",
    "test": "*PR.wdl"
}

class TestFailure(Exception):
    """Error for a test that fails its portability test."""
    pass

class TestDefinitionError(Exception):
    """Error for malformed test directory."""
    pass

class TestEnvironmentError(Exception):
    """Error for malformed test environment."""
    pass

def gather_test_inputs(test_dir):
    """Walk through the test directory, finding files and dependencies needed
    to submit a portability test.

    Raise a TestDefinitionError if anything goes wrong.
    """

    errors = []

    # Test WDL
    wdl_glob = os.path.join(test_dir, TEST_DIR_LAYOUT["test"])
    wdl_paths = glob.glob(wdl_glob)
    if not wdl_paths:
        errors.append("Test definition WDL not found.")
    if len(wdl_paths) > 1:
        errors.append("Multiple candidate test WDLs found: {}".format(wdl_paths))

    workflow_attachment = []
    try:
        test_wdl_name = os.path.basename(wdl_paths[0])
        test_wdl_string = open(wdl_paths[0]).read()
        workflow_attachment.append((test_wdl_name, test_wdl_string))
    except IOError:
        test_wdl_name, test_wdl_string = None, None
        errors.append("Test WDL {} could not be read".format(wdl_paths[0]))

    # Test inputs
    try:
        inputs_json_path = os.path.join(test_dir, TEST_DIR_LAYOUT["inputs"])
        inputs_json_string = open(inputs_json_path).read()
    except IOError:
        inputs_json_string = None
        errors.append("Inputs JSON {} could not be read".format(inputs_json_path))

    # Dependencies

    # First try to load the dependency JSON itself
    try:
        dependencies_json_path = os.path.join(test_dir, TEST_DIR_LAYOUT["dependencies"])
        dependencies_json_string = open(dependencies_json_path).read()
    except IOError:
        dependencies_json_string = None
        errors.append("Dependencies JSON {} could not be read".format(dependencies_json_path))

    # Now iterate over the dependencies and load the files.
    dependencies_dict = json.loads(dependencies_json_string)

    for key, value in dependencies_dict.items():
        try:
            workflow_attachment.append((
                key,
                open(value).read()
            ))
        except IOError:
            errors.append("Could not read dependency {}".format(value))

    if errors:
        for error in errors:
            sys.stderr.write(error + "\n")
        raise TestDefinitionError(errors)

    return {
        "entry_point_wdl": test_wdl_name,
        "workflow_params": inputs_json_string,
        "workflow_attachment": workflow_attachment
        }

# These are the environment variables we're expecting to be able to submit a
# portability test.
ENV_VARIABLES = {
    "portability_service_url": os.environ.get("PORTABILITY_SERVICE_URL"),
    "portability_service_headers": os.environ.get("PORTABILITY_SERVICE_HEADERS")
}

def verify_environment_variables():
    """Check that required environment variables are defined."""

    errors = []

    for key, value in ENV_VARIABLES.items():
        if not value:
            errors.append("Environment variable {} is undefined".format(key.upper()))
    if errors:
        for error in errors:
            sys.stderr.write(error + "\n")
        raise TestEnvironmentError(errors)

def submit_portability_test(test_inputs):
    """Submit the portability to the service and return the test's id."""

    service_headers = json.loads(ENV_VARIABLES["portability_service_headers"])
    test_endpoint = ENV_VARIABLES["portability_service_url"] + "/portability_tests"
    response = requests.post(
        test_endpoint,
        headers=service_headers,
        json=test_inputs)

    print("Portability service response:\n{}".format(
        json.dumps(json.loads(response.text), indent=4)), flush=True)

    test_id = json.loads(response.text)["test_id"]

    return test_id

def print_test_states(test_states):
    """Print the states of tests."""
    print("{:%Y-%m-%dT%H:%M:%SZ}".format(datetime.datetime.utcnow()), flush=True)
    for test_name, test_state in test_states.items():
        print("TEST:{} STATE:{}".format(test_name, test_state), flush=True)

def monitor_tests(test_ids):
    """Check the status of tests until they all either fail or succeed. If any fail, raise
    TestFailure, and if they all succeed, return successfully.
    """

    service_headers = json.loads(ENV_VARIABLES["portability_service_headers"])

    terminal_tests = {}
    while True:
        time.sleep(120)

        test_states = {}

        for test_name, test_id in test_ids.items():

            if test_name in terminal_tests:
                test_states[test_name] = terminal_tests[test_name]
                continue

            status_endpoint = ENV_VARIABLES["portability_service_url"]  + "/portability_tests/" + \
                test_id + "/status"

            response = requests.get(status_endpoint, headers=service_headers)
            test_state = json.loads(response.text)["state"]
            test_states[test_name] = test_state

        print_test_states(test_states)

        for test_name, test_state in test_states.items():
            if test_state in ("Failed", "Succeeded"):
                terminal_tests[test_name] = test_state

        if len(terminal_tests) == len(test_ids):
            break

    for test_name, test_id in test_ids.items():

        if terminal_tests[test_name] != "Succeeded":
            info_endpoint = ENV_VARIABLES["portability_service_url"]  + "/portability_tests/" + \
                test_id
            response = requests.get(info_endpoint, headers=service_headers)
            print(response.json())

    if not all(k == "Succeeded" for k in terminal_tests.values()):
        raise TestFailure()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--test-directories", required=True, nargs='+')
    args = parser.parse_args()

    verify_environment_variables()

    tests_inputs = {}
    for test_directory in args.test_directories:
        test_inputs = gather_test_inputs(test_directory)
        tests_inputs[test_directory] = test_inputs

    test_ids = {}
    for test_dir, test_inputs in tests_inputs.items():
        test_id = submit_portability_test(test_inputs)
        test_ids[test_dir] = test_id

    monitor_tests(test_ids)

if __name__ == "__main__":
    main()
