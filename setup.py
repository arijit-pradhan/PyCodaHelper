import setuptools

setuptools.setup(
    include_package_data=True,
    name='PyCodaHelper',
    version='0.0.1',
    description='CoDA Helper Python Module',
    #url='...',
    author='-',
    author_email='...@wsu.edu',
    packages=setuptools.find_packages(),
    install_requires=['numpy', 'vtk'],
    long_description='Helper package, can be used with CoDA',
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
         "Operating System :: OS Independent",
    ],
)