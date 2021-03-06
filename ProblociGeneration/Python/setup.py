# LICENSE
# Copyright (c) 2018-2021 Genome Research Ltd.
# Author: CASM-IT <cgphelp@sanger.ac.uk>
#
#
# This file is part of battenberg.
#
# problocifileparser is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#    1. The usage of a range of years within a copyright statement contained within
#    this distribution should be interpreted as being equivalent to a list of years
#    including the first and last year specified and all consecutive years between
#    them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
#    2009, 2011-2012’ should be interpreted as being identical to a statement that
#    reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
#    statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
#    identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
#    2009, 2010, 2011, 2012’."
#
#

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    "description": "Parse prob loci initial output to provide final prob loci file",
    "author": "David Jones",
    "url": "https://gitlab.internal.sanger.ac.uk/CancerIT/battenberg_generate_probloci/",
    "download_url": "https://gitlab.internal.sanger.ac.uk/CancerIT/battenberg_generate_probloci/1.0.0.tar.gz",
    "author_email": "cgpit@sanger.ac.uk",
    "version": "1.0.0",
    "install_requires": ["pysam", 
                          "numpy", 
                          "matplotlib"],
    "packages": ["problocifileparser"],
    "scripts": ["parse_prob_loci_stats.py"],
    "name": "problocifileparser",
}

setup(**config)
