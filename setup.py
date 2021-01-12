import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cwipp",
    version="0.0.1",
    author="Nadya Moisseeva",
    author_email="nmoisseeva@eoas.ubc.ca",
    description="Plume-rise parameterization package for wildfire smoke",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nmoisseeva/cwipp/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
