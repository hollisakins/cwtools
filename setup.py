from setuptools import find_packages, setup

setup(
    name='cwtools',
    packages=find_packages(include=['cwtools']),
    version='0.1.0',
    description='Tools for COSMOS-Web',
    author='Hollis Akins',
    license='MIT',
    install_requires=['astropy','numpy','matplotlib','requests','aiohttp','rich','pwinput','tqdm','requests'],
)
