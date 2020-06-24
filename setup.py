# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='percival',
    version='0.1.11',
    description='input generator for quantum chemical calculation',
    long_description=readme,
    author='Yusuke Kawashima',
    author_email='kawashima-y@phs.osaka-u.ac.jp',
    install_requires=[],
    include_package_data=True,
    url='https://github.com/kyusque/percival',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    test_suite='tests'
)

