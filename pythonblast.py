"""
**pythonblast.py
** Copyright (C) 2018  Jose Sergio Hleap

Python wrapper for multithreaded blast one sequence at a time

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: jose.hleaplozano@mcgill.ca
"""
from io import StringIO
import pandas as pd
from multiprocessing import Pool
from joblib import Parallel, delayed
from subprocess import Popen, PIPE
import optparse
from itertools import zip_longest
import gzip
import csv
import os
