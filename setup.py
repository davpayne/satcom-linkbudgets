import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='satcom_linkbudgets',
      version='0.1',
      description='Satellite Link Budgets',
      author="David Payne",
    author_email="dpayne162@gmail.com",
    description="A package for calculating satellite link budgets, optical and RF",
    url="https://github.com/davpayne/satcom-linkbudgets",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = ['numpy'],
    python_requires='>=3.6',
      packages=['link_tools'],
      zip_safe=False)
