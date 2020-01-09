from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
    long_description = readme_file.read()

setup(
    name="Heat-Content",
    version="1.0.0",
    description="Python packge to calculate Upper Ocean Heat Content from ROMS model output",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Patrick Daniel",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License"
    ],
    packages=find_packages(),
    python_requires=">3.5"
)