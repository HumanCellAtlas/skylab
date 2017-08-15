import setuptools

setuptools.setup(
    name="martian-cli",
    version="0.0.1",
    description="An alternate command line interface to 10x Genomics Martian workflows.",
    url="https://github.com/HumanCellAtlas/skylab/",
    author="Marcus Kinsella",
    author_email="mkinsella@chanzuckerberg.com",
    packages=["martian_cli"],
    install_requires=["pyparsing"],
    entry_points={
        'console_scripts': ["martian=martian_cli.cli:main"]
        }
)
