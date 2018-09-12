import sys
import subprocess as sp
import os

bams_dir = sys.argv[1]
r5_index_file = sys.argv[2]
new_bams_dir = sys.argv[3]

bam_files = str.strip(sp.check_output("ls " + bams_dir + "/*.unique*bam", shell = True)).split("\n")

r5_indexes = set()
with open(r5_index_file) as f:
    for line in f:
        r5_indexes.add(str.strip(line))

bam_files_move = []
for fi in bam_files:
    fields = str.strip(fi).split("_")
    r5 = fields[-1].split(".")[0]
    if r5 in r5_indexes:
        bam_files_move.append(fi)

os.mkdir(new_bams_dir)
for fi in bam_files_move:
    mv_cmd = "mv " + fi + " " + new_bams_dir
    sp.call(mv_cmd, shell = True)

