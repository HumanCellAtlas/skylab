## Running the IPython notebook

To run 10x_cost_analysis, first ensure you have dependencies properly installed:

`pip3 install -r requirements.txt`

Then, run the notebook as follows, passing in the required metadata parameter

METADATA=your_filename.json jupyter nbconvert --execute --inplace \
--to notebook 10x_cost_analysis.ipynb
