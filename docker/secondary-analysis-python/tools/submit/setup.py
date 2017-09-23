from setuptools import setup

setup(name='submit',
      version='1.0.0.dev1',
      description='submit',
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
              'submit=submit.submit:main',
              'analysis-json=submit.analysis_json:main',
              'stage=submit.stage:run',
              'submission-urn=submit.submission_urn:main',
              'confirm=submit.confirm:main'
          ]
      },
      include_package_data=True
)
