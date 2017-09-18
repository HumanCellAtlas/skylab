from setuptools import setup

setup(name='submit',
      version='1.0.0.dev1',
      description='submit',
      url='http://github.com/HumanCellAtlas/skylab',
      author='Dave Shiga',
      author_email='dshiga@broadinstitute.org',
      license='MIT',
      packages=['submit', 'tests'],
#      package_dir={'': 'src'},
      entry_points = {
          "console_scripts": ['hca-submit=submit.submit:main']
      },
      include_package_data=True
)
