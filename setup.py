#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from setuptools import setup, find_packages
import pathlib
here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")
with open("README.md", "r", encoding='utf-8') as fh:
    long_des = fh.read()

required = ['pandas', 'numpy', 'biopython', 'matplotlib', 'scipy', 'seaborn', 'plotnine', 'ete3', 'scikit-learn']

setup(
    name="quota_anchor",
    version="1.0.0",
    author="XiaoDong Li",
    author_email="xiaodongli2405@gmail.com",
    description="Conduct strand and WGD aware syntenic identification",
    license="MIT License",
    long_description=long_des,
    long_description_content_type="text/markdown",
    url="https://github.com/baoxingsong/quota_Anchor",
    packages=['quota_anchor', 'quota_anchor.lib', 'quota_anchor.kspeaks', 'quota_anchor.config_file', 'quota_anchor.plots', 'scripts'],
    package_data={
        'quota_anchor': ['config_file/*', 'plots/*'],
    },
    include_package_data=True,
    python_requires=">=3, <4",
    classifiers=[
        "Development Status :: 3 - RC",
        "Intended Audience :: Comparative genomics/Evolution researcher",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
    ],
    entry_points={
        'console_scripts': [
            'quota_Anchor = quota_anchor.main:main',
        ]
    },
    zip_safe=True,
    install_requires=required
)
