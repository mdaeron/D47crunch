import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="D47crunch-mdaeron", # Replace with your own username
    version="0.0.1",
    author="Mathieu Daëron",
    author_email="daeron@lsce.ipsl.fr",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mdaeron/D47crunch",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Modified BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)