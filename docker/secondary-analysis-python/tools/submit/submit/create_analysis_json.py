#!/usr/bin/env python

import requests
import json
import argparse
from collections import OrderedDict

def create_analysis(analysis_id, metadata_file, input_bundles_string, reference_bundle,
        run_type, method, schema_version, inputs_file, outputs_file, format_map):
    print('Creating analysis.json for {}'.format(analysis_id))

    with open(metadata_file) as f:
        metadata = json.load(f)
        start, end = get_start_end(metadata)
        tasks = get_tasks(metadata)

    inputs = create_inputs(inputs_file)
    outputs = create_outputs(outputs_file, format_map)

    input_bundles = get_input_bundles(input_bundles_string)

    analysis = {}
    analysis['analysis_id'] = analysis_id
    analysis['analysis_run_type'] = 'run'
    analysis['reference_bundle'] = reference_bundle
    analysis['computational_method'] = method
    analysis['input_bundles'] = input_bundles
    analysis['timestamp_start_utc'] = start
    analysis['timestamp_stop_utc'] = end
    analysis['metadata_schema'] = schema_version
    analysis['tasks'] = tasks
    analysis['inputs'] = inputs
    analysis['outputs'] = outputs

    return analysis

def create_inputs(inputs_file):
    inputs = []
    with open(inputs_file) as f:
        f.readline() # skip header
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            value = parts[1]
            input = {
                'name': name,
                'value': value
            }
            if value.startswith('gs://'):
                # This is a placeholder for now, since the analysis json schema requires it.
                # In future we will either properly calculate a checksum for the file
                # or remove it from the schema.
                input['checksum'] = 'd0f7d08f1980f7980f'
            inputs.append(input)
    return inputs

def create_outputs(outputs_file, format_map):
    with open(format_map) as f:
        extension_to_format = json.load(f)

    outputs = []
    with open(outputs_file) as f:
        for line in f:
            path = line.strip()
            d = {
              'file_path': path,
              'name': path.split('/')[-1],
              'format': get_format(path, extension_to_format)
            }
            outputs.append(d)
    return outputs

def get_format(path, extension_to_format):
    for ext in extension_to_format:
        if path.endswith(ext):
            format = extension_to_format[ext]
            print(format)
            return format
    print('Warning: no known format matches file {}'.format(path))
    return 'unknown'

def get_input_bundles(input_bundles_string):
    input_bundles = input_bundles_string.split(',')
    print(input_bundles)
    return input_bundles

def get_start_end(metadata):
    start = metadata['start']
    end = metadata['end']
    print(start, end)
    return start, end

def get_tasks(metadata):
    calls = metadata['calls']
    out_tasks = []
    for long_task_name in calls:
        out_task = {}
        task_name = long_task_name.split('.')[-1]
        task = calls[long_task_name][0]
        out_task['name'] = task_name
        runtime = task['runtimeAttributes']
        out_task['cpus'] = int(runtime['cpu'])
        out_task['memory'] = runtime['memory']
        out_task['disk_size'] = runtime['disks']
        out_task['docker_image'] = runtime['docker']
        out_task['zone'] = runtime['zones']
        out_task['start_time'] = task['start']
        out_task['stop_time'] = task['end']
        out_task['log_out'] = task['stdout']
        out_task['log_err'] = task['stderr']
        out_tasks.append(out_task)
    sorted_out_tasks = sorted(out_tasks, key=lambda k: k['name'])
    return sorted_out_tasks

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-analysis_id', required=True)
    parser.add_argument('-metadata_json', required=True)
    parser.add_argument('-input_bundles', required=True)
    parser.add_argument('-reference_bundle', required=True)
    parser.add_argument('-run_type', required=True)
    parser.add_argument('-method', required=True)
    parser.add_argument('-schema_version', required=True)
    parser.add_argument('-inputs_file', required=True)
    parser.add_argument('-outputs_file', required=True)
    parser.add_argument('-format_map', required=True)
    args = parser.parse_args()
    analysis = create_analysis(args.analysis_id, args.metadata_json, args.input_bundles,
        args.reference_bundle, args.run_type, args.method, args.schema_version,
        args.inputs_file, args.outputs_file, args.format_map)
    with open('analysis.json', 'w') as f:
        json.dump(analysis, f, indent=2, sort_keys=True)

if __name__ == '__main__':
    main()
