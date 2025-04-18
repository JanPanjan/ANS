import sys
from time import time

if len(sys.argv) != 2:
	raise Exception("Usage: python3 lines.py <fastq-file>")

fname = sys.argv[1]

print("Processing file...")
start_time = time()

with open(fname, "r") as file:
    data = file.readlines()
    n = len(data) / 4

end_time = time()
final_time = end_time - start_time

print(f"number of sequences in file: {n}")
print(f"time: {final_time} sec")
