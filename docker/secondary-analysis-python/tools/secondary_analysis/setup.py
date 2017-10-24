from setuptools import setup

setup(name='secondary_analysis',
      version='1.0.0.dev1',
      description='Utilities for retrieving files from the HCA data storage service and submitting an analysis bundle to HCA-DCP',
      url='http://github.com/HumanCellAtlas/skylab',
      author='Dave Shiga',
      author_email='dshiga@broadinstitute.org',
      license='BSD 3-clause "New" or "Revised" License',
      packages=['secondary_analysis', 'tests'],
      install_requires=[
          'requests',
          'boto3'
      ],
      entry_points={
          "console_scripts": [
              'create-analysis-json=secondary_analysis.create_analysis_json:main',
              'create-envelope=secondary_analysis.create_envelope:main',
              'get-staging-urn=secondary_analysis.get_staging_urn:main',
              'stage=secondary_analysis.stage:run',
              'confirm-submission=secondary_analysis.confirm_submission:main'
          ]
      },
      include_package_data=True
      )
