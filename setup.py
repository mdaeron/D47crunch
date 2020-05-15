import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="D47crunch",
    version="0.3.2",
    author="Mathieu Daëron",
    author_email="daeron@lsce.ipsl.fr",
    description="Standardization of Δ47 clumped-isotope measurements",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mdaeron/D47crunch",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
			'lmfit >= 1.0.1',
			'matplotlib >= 3.1.3',
			'numpy >= 1.18.1',
			'scipy >= 1.4.1',
    	]
)