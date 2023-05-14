# standard files
import os
from glob import glob
import sys

for directory in glob(f"{sys.argv[1]}/*"):
    runs=glob(f"{directory}/*")
    for run in runs:
        if not len(glob(f"{run}/final*")):
            for file in os.listdir(run):
                os.remove(f"{run}/{file}")
