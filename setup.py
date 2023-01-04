import setuptools
import versioneer

with open("README.rst", "r") as fh:
    long_description = fh.read()
with open("requirements.txt", "r") as fh:
    requirements = [line.strip() for line in fh]

setuptools.setup(
    name="reco",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Martin Wegner",
    author_email="martinwegner1983@gmail.com",
    description="A Python library to find gRNA read counts in fastq files.",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requirements=requirements,
)
