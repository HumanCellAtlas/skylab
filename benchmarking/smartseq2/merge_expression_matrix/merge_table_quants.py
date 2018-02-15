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
run_name=sys.argv[2]
output_name=sys.argv[4]
value_name=sys.argv[3]
##meta_url
metadata_url="https://cromwell.mint-dev.broadinstitute.org/api/workflows/v1/"+uuid+"/metadata?expandSubWorkflows=false"
r=requests.get(metadata_url,auth=(logins['cromwell_username'],logins['cromwell_password']))
data=r.json()
## load output files
files=data['outputs'][run_name]
## parsing
merged=pd.DataFrame()
for kk in range(0,len(files)):
    print(kk)
    fc1=files[kk]
    ##samplename=data[kk]['inputs']['sample_name']
    fc1=fc1.replace('gs://broad-dsde-mint-dev-cromwell-execution/','')
    blob1=bucket.get_blob(fc1)
    bname1=basename(fc1)
    print(bname1)
    sample_name=bname1.split('.')[0]
    with open(bname1,'wb') as file_obj:
        blob1.download_to_file(file_obj)
    if "unq_genes_counts" in run_name or "mult_genes_counts" in run_name:
        if kk==0:
            merged=pd.read_csv(bname1,skiprows=2,sep='\t',usecols=['gene_id','length',sample_name],names=['gene_id','contigs','start','end','strands','length',sample_name])
            if value_name == 'tpm': ## convert from cnt to tpm
                sample=list(merged[sample_name].values)
                lengths=list(merged['length'].values)
                rpk=[float(s) / float(l) for s,l in zip(sample, lengths)]
                RPKM=sum(rpk)/1000000
                TPM=[r/RPKM for r in rpk]
                merged[sample_name]=TPM
        else:
            cnt=pd.read_csv(bname1,skiprows=2,sep='\t',usecols=['gene_id','length',sample_name],names=['gene_id','contigs','start','end','strands','length',sample_name])
            if value_name == "tpm":
                sample=list(cnt[sample_name].values)
                lengths=list(cnt['length'].values)
                rpk=[float(s) / float(l) for s,l in zip(sample, lengths)]
                RPKM=sum(rpk)/1e6
                TPM=[r/RPKM for r in rpk]
                cnt[sample_name]=TPM
            merged=pd.merge(left=merged,right=cnt[['gene_id',sample_name]],left_on='gene_id',right_on='gene_id')
    elif 'rsem' in run_name:
        if kk==0:
            merged=pd.read_csv(bname1,delimiter='\t',skiprows=1,sep='\t',names=['gene_id','trans_id', 'length', 'eff_length','est_counts','tpm','fpkm','pme_count','pme_sd','pme_tpm','pme_fpkm'],usecols=['gene_id','length',value_name])
            merged.columns=['gene_id','length',sample_name]
        else:
            cnt=pd.read_csv(bname1,delimiter='\t',skiprows=1,sep='\t',names=['gene_id','trans_id', 'length', 'eff_length','est_counts','tpm','fpkm','pme_count','pme_sd','pme_tpm','pme_fpkm'],usecols=['gene_id',value_name])
            cnt.columns=['gene_id',sample_name]
            merged=pd.merge(left=merged,right=cnt,left_on='gene_id',right_on='gene_id')
merged=merged.round(3)
merged.to_csv(output_name,index=False)

