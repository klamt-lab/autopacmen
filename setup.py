import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="autopacmen-Paulocracy",
    version="0.6.0",
    author="Paulocracy",
    author_email="bekiaris@mpi-magdeburg.mpg.de",
    description="The AutoPACMEN package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Paulocracy/autopacmen",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    install_requires=[
        "biopython",
        "cobra",
        "click",
        "openpyxl",
        "pebble",
        "requests",
        "xlsxwriter",
    ],
)
