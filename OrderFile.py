import sys
import numpy as np

filename_in = str(sys.argv[1])
if len(sys.argv) == 3:
    filename_out = str(sys.argv[2])
else:
    filename_out = filename_in


with open(filename_in, "r") as f:
    sorted_file = sorted(f)

with open(filename_out, "w+") as f:
    for line in sorted_file:
        f.write(line)
