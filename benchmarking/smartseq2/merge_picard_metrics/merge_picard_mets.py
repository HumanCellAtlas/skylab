from crimson import picard
import pandas as pd
import numpy as np
from google.cloud import storage
import json
from os.path import basename
import sys
import requests

## set up auth
client = storage.Client()
bucket=client.get_bucket('broad-dsde-mint-dev-cromwell-execution')
## load cromwell credential
logins=json.load(open('/usr/secrets/broad-dsde-mint-dev-cromwell.json'))

## input json has uuid, run name and value to parse.
uuid=sys.argv[1]
met_name=sys.argv[2]
output_name=sys.argv[3]
##meta_url
metadata_url="https://cromwell.mint-dev.broadinstitute.org/api/workflows/v1/"+uuid+"/metadata?expandSubWorkflows=false"
print(metadata_url)
r=requests.get(metadata_url,auth=(logins['cromwell_username'],logins['cromwell_password']))
data=r.json()
## load output files
files=data['outputs'][met_name]
## parsing
mets={}
for kk in range(0,len(files)):
    fc1=files[kk]
    fc1=fc1.replace('gs://broad-dsde-mint-dev-cromwell-execution/','')
    blob1=bucket.get_blob(fc1)
    bname1=basename(fc1)
    sample_name=bname1.split('.')[0]
    with open(bname1,'wb') as file_obj:
        blob1.download_to_file(file_obj)
    parsed = picard.parse(bname1)
    class_name = parsed['metrics']['class']
    if class_name == "picard.analysis.AlignmentSummaryMetrics":
        ## only parse out pair reads 
        met = parsed['metrics']['contents'][2]
    elif class_name == "picard.analysis.InsertSizeMetrics":
        if len(parsed['metrics']['contents']) <21:
            met = parsed['metrics']['contents'][0]
        else:
            met = parsed['metrics']['contents']
    else:
        met = parsed['metrics']['contents']
    mets[sample_name] = met
tab = pd.DataFrame.from_dict(mets)
tab.to_csv(output_name)

