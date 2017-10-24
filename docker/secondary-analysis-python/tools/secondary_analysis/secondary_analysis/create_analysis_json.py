#!/usr/bin/env python

import json
import argparse

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

    analysis = {
        'analysis_id': analysis_id,
        'analysis_run_type': run_type,
        'reference_bundle': reference_bundle,
        'computational_method': method,
        'input_bundles': input_bundles,
        'timestamp_start_utc': start,
        'timestamp_stop_utc': end,
        'metadata_schema': schema_version,
        'tasks': tasks,
        'inputs': inputs,
        'outputs': outputs
    }

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
        task_name = long_task_name.split('.')[-1]
        task = calls[long_task_name][0]
        runtime = task['runtimeAttributes']

        out_task = {
            'name': task_name,
            'cpus': int(runtime['cpu']),
            'memory': runtime['memory'],
            'disk_size': runtime['disks'],
            'docker_image': runtime['docker'],
            'zone': runtime['zones'],
            'start_time': task['start'],
            'stop_time': task['end'],
            'log_out': task['stdout'],
            'log_err': task['stderr']
        }
        out_tasks.append(out_task)
    sorted_out_tasks = sorted(out_tasks, key=lambda k: k['name'])
    return sorted_out_tasks

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysis_id', required=True, help='Cromwell workflow id')
    parser.add_argument('--metadata_json', required=True, help='JSON obtained from calling Cromwell metadata endpoint for this workflow id')
    parser.add_argument('--input_bundles', required=True, help='The DSS bundles used as inputs for the workflow')
    parser.add_argument('--reference_bundle', required=True, help='To refer to the DSS resource bundle used for this workflow, once such things exist')
    parser.add_argument('--run_type', required=True, help='Should always be "run" for now, may be "copy-forward" in some cases in future')
    parser.add_argument('--method', required=True, help='Supposed to be method store url, for now can be url for wdl in skylab')
    parser.add_argument('--schema_version', required=True, help='The metadata schema version that this analysis.json conforms to')
    parser.add_argument('--inputs_file', required=True, help='Path to tsv containing info about inputs')
    parser.add_argument('--outputs_file', required=True, help='Path to json file containing info about outputs')
    parser.add_argument('--format_map', required=True, help='JSON file providing map of file extensions to formats')
    args = parser.parse_args()
    analysis = create_analysis(args.analysis_id, args.metadata_json, args.input_bundles,
        args.reference_bundle, args.run_type, args.method, args.schema_version,
        args.inputs_file, args.outputs_file, args.format_map)
    with open('analysis.json', 'w') as f:
        json.dump(analysis, f, indent=2, sort_keys=True)

if __name__ == '__main__':
    main()
