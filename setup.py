#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = ['pytest>=3', ]

setup(
    author="Anastasiia Kornilova",
    author_email='anastasiia.kornilova@skoltech.ru',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description=" Metrics to evaluate map and trajectory consistency via estimating quality of map aggregated from point clouds",
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description=readme + '\n\n' + history,
    long_description_content_type='text/markdown'
    include_package_data=True,
    keywords='map_metrics',
    name='map_metrics',
    packages=find_packages(include=['map_metrics', 'map_metrics.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/anastasiia-kornilova/map_metrics',
    version='0.0.1',
    zip_safe=False,
)
