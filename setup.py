import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GPDynamics",
    version="0.0.1",
    author="Marcio Gameiro",
    author_email="marciogameiro@gmail.com",
    description="Gaussian Process Dynamics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/marciogameiro/GPDynamics",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.0',
)
