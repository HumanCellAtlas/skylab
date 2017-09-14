# Overview

This directory contains a wrapper WDL used by the HCA Data Coordination Platform's Secondary Analysis Service to run the ss2_single_sample WDL. It is provided here mainly for informational purposes. Most people who want to run the SmartSeq2 pipeline should just run the ss2_single_sample.wdl directly.

# Files

* wrapper_ss2_single_sample.wdl
The wrapper wdl, which calls prepare_and_analyze_ss2_single_sample.wdl and then submit_ss2_single_sample.wdl

* prepare_and_analyze_ss2_single_sample.wdl
Given a bundle uuid and version as inputs, queries the DCP Storage Service and obtain the actual fastq gs URLs needed to run ss2_single_sample.wdl, which it then runs.

* submit_ss2_single_sample.wdl
Submits the results from ss2_single_sample.wdl back to the DCP Ingest Service.

# Do I need this?

Only the Secondary Analysis Service should be trying to submit SmartSeq2 outputs to DCP Ingest. As mentioned above, if you want to run the SmartSeq2 pipeline you should probably run the ss2_single_sample.wdl directly. The prepare_and_analyze_ss2_single_sample.wdl might be useful to you if you want to translate bundle uuids and versions to fastq gs URLs.
