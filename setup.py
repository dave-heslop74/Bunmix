import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bunmix",
    version="1.0.0",
    author="D. Heslop",
    author_email="notme@gmail.com",
    description="The Bunmix package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Bunmix",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
