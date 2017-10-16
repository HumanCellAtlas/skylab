This directory contains WDLs used by the Secondary Analysis Service to interface with the Human Cell Atlas Data Storage Service and Ingestion Service.

# submit.wdl

This WDL is used to create analysis bundles. It talks to the Ingestion Service, which in turn creates bundles in the Data Storage Service.

For more information, see comments in submit.wdl and the README.rst in the "submit" Python package in this repo, located at:
/docker/secondary-analysis-python/tools/submit
