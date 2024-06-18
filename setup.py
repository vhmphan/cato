from setuptools import setup, find_packages

setup(
    name='cato',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy'
    ],
    description='A package containing a function to cosmic-ray ionization cross-section',
    author='Vo Hong Minh Phan',
    author_email='vhmphan@obspm.fr',
    url='https://github.com/vhmphan/cato',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)