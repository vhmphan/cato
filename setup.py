from setuptools import setup, find_packages

setup(
    name='cato',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy'
    ],
    description='A package containing a function to calculate electron capture cross-section',
    author='Your Name',
    author_email='your.email@example.com',
    url='https://github.com/yourusername/mypackage',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)