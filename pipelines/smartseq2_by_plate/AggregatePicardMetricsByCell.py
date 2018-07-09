from crimson import picard
import json
import argparse

def AggregatePicardMetricByCell(files, output_name):
    """
    piepline output picard QC metrics at sinle cell/sample level.
    This functin is called to merge/aggregate QC metrics by metrics type and then merge multiple QC measurement 
    into single matrix file. In this file, column is sample/cell and row is QC metrics
    :param files: metric files from pipeline outputs
    :param met_name: metrics name with workflow name and subworkflow name as prefix. such as 'run_pipelines.RunStarPipeline.alignment_summary_metrics'
    """
    # initial output
    print(files)
    filenames = files.split(',')
    for kk in range(0, len(filenames)):
        file_name = filenames[kk]
        parsed = picard.parse(file_name)
        class_name = parsed['metrics']['class']
        # Aignment metrics return multiple lines, but only output PAIRED-READS/third line
        if class_name == "picard.analysis.AlignmentSummaryMetrics":
            ## only parse out pair reads
            met = parsed['metrics']['contents'][2]
        # sometimes(very rare), insertion metrics also return multiple lines results to include TANDEM repeats. but we only output the first line.
        elif class_name == "picard.analysis.InsertSizeMetrics":
            # if the elemnet counts is less than 21, it means insertion metrics returns multiple line results.
            if len(parsed['metrics']['contents']) < 21:
                met = parsed['metrics']['contents'][0]
            else:
                met = parsed['metrics']['contents']
        else:
            # other metrics(so far) only return one line results.
            met = parsed['metrics']['contents']
        if kk == 0 :
            mets=met.copy()
        else:
            mets.update(met)
    f=open(output_name+'.picard.core.json','w')
    f.write(json.dumps({output_name:mets},indent=2))
    f.close()
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--file_names",
        dest="file_names",
        required=True,
        help=
        "a list of files to be parsed out."
    )
    parser.add_argument(
        "-o",
        "--output_name",
        dest="output_name",
        required=True,
        help="The output file name")
    args = parser.parse_args()
    AggregatePicardMetricByCell(args.file_names, args.output_name)

if __name__ == "__main__":
    main()
