import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="formulate", 
    packages=setuptools.find_packages(include=['formulate']),
    version="1.1",
    author="Luc E.Brunet",
    author_email="luc.brunet@rd-mediation.fr",
    description="Library for formulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://www.rd-mediation.fr",
    #packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Proprieraty",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
) 