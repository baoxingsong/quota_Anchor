#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from setuptools import setup, find_packages
import pathlib
here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")
with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

required = ['pandas', 'numpy', 'biopython', 'matplotlib', 'scipy', 'seaborn', 'plotnine']

setup(
    name="quota_Anchor",
    version="0.0.1_alpha",
    author="XiaoDong Li",
    author_email="xiaodongli2405@gmail.com",
    description="Conduct strand and WGD aware syntenic identification",
    license="MIT License",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/baoxingsong/quota_Anchor",
    packages=find_packages(),
    package_data={'quota_Anchor': ['config_file/*']},
    python_requires=">=3, <4",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Comparative genomics/Evolution researcher",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
    ],
    entry_points={
        'console_scripts': [
            'quota_Anchor = quota_Anchor.main:main',
        ]
    },
    zip_safe=True,
    install_requires=required
)