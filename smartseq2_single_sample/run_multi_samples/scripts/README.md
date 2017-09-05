# Running Time and Cost Calculation
This script is original post in [dsde-pipeline](https://github.com/broadinstitute/dsde-pipelines/blob/develop/scripts/calculate_cost.py), which is used to calculate running time and cost of workflow in google cloud. 
We modified this python script in order to calculate cost of a subworkflow.
# Requirement
- Python3
- Required modules
 --------
 ```python
 import argparse
import json
import math
import urllib.request, urllib.error, urllib.parse
import dateutil.parser
from io import StringIO
import gzip
import requests
```
# Usage and Basic Options
```bash
usage: calculate_cost_workflow.py [-h] -i INPUT_JSON_FILE

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_JSON_FILE, --input-json INPUT_JSON_FILE
                        input metadata information
```
# Inputs
```json
{
  "username": "login-name of cromwell",
  "password": "login password",
  "master_metadata":"Json metadata output of workflow",
  "workflow": "subworkflow name or aliase",
  "run_identifier":"the id to identify each run in subworkflow,such as unique sample id, cell id ",
  "ignore_preempted": 0,
  "only_total_cost": 0
}

```
