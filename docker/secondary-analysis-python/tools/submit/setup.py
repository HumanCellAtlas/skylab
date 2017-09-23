from setuptools import setup

setup(name='submit',
      version='1.0.0.dev1',
      description='Utility for submitting an analysis bundle to HCA-DCP',
      url='http://github.com/HumanCellAtlas/skylab',
      author='Dave Shiga',
      author_email='dshiga@broadinstitute.org',
      license='MIT',
      packages=['submit', 'tests'],
      install_requires=[
          'requests',
      ],
      entry_points = {
          "console_scripts": [
              'create-analysis-json=submit.create_analysis_json:main',
              'create-envelope=submit.create_envelope:main',
              'get-staging-urn=submit.get_staging_urn:main',
              'stage=submit.stage:run',
              'confirm-submission=submit.confirm_submission:main'
          ]
      },
      include_package_data=True
)
