from sample_histoandmeastd import *
from config import *
from mulit_2lumpy import *
from BayesianOptimizer import *
import os
import subprocess
import tempfile
import config
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
from tqdm import tqdm
import pandas as pd
import re
from multiprocessing import cpu_count, Pool
    


sample_histoandmeastdmain()
mulit_2lumpymain()
BayesianOptimizermain()