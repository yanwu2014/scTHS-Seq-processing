import pysam as ps
import sys
from collections import Counter
from itertools import product
import subprocess as sp
import multiprocessing as mp
import os
import time
from functools import partial


def main():
    read1_fastq_file = sys.argv[1]
    i7_index_fastq_file = sys.argv[2]
    i5_index_fastq_file = sys.argv[3]
    out_prefix = sys.argv[4]
    bwa_index_path = sys.argv[5]
    n_cores = int(sys.argv[6])
    
    ## Time the script
    start = time.time()

    ## Hardcoded parameters
    i7_index_list = "i7_scTHS"
    i5_index_list = "i5_scTHS"
    r5_index_list = "r5_scTHS"
    n_mismatch = 1
    min_unmapped_reads = 500


    ## Split out index files and generate config 
    print "Splitting index files and generating config file..."
    split_index(i7_index_fastq_file, out_prefix + ".i7", 0, 8)
    split_index(i5_index_fastq_file, out_prefix + ".i5", 0, 8)
    split_index(i5_index_fastq_file, out_prefix + ".r5", 23, 29)
    
    ## Build config file for deindexer
    build_config(i7_index_list, i5_index_list, r5_index_list, "scTHS")
    
    ## Run deindexer
    print "Running fastq deindexing..."
    run_deindexer(out_prefix, read1_fastq_file, i7_index_list, i5_index_list, r5_index_list)
    
    ## Merge fastq files and stash cell index in read name
    print "Merging fastq files for alignment..."
    cell_read_counts = merge_fastq(out_prefix + "_deindexed_fastq", out_prefix)
    
    ## Align merged fastq files
    print "Running alignment with BWA..."
    run_alignment(out_prefix, bwa_index_path, n_cores)

    ## Deindex bams using stashed cell indexes
    print "Deindexing aligned reads using stashed cell indexes..."
    valid_cells = [k for k,v in cell_read_counts.items() if v >= min_unmapped_reads]
    print "Cell Barcodes Passing Filter:\t" + str(len(valid_cells))
    unmerge_bams(out_prefix, valid_cells, n_cores)
    
    ## Remove PCR duplicates
    dedup_bams(out_prefix + "_deindexed_bam", n_cores)
    
    ## Cleanup some temproary files
    sp.call("rm -r " + out_prefix + "_deindexed_fastq/", shell = True)
    sp.call("rm *.fastq", shell = True)

    print "Time elapsed:\t" + str(time.time() - start)



def count_reads(bams_dir):
    deduped_files = str.strip(sp.check_output("ls " + bams_dir + "/*.unique.bam", shell = True)).split("\n")
    unique_reads = {}
    for fi in deduped_files:
        unique_reads[fi] = ps.AlignmentFile(fi, "rb").mapped

    return unique_reads


def _dedup_index_bam_(filename):
    rmdup_cmd = "samtools rmdup -s " + filename + " " + filename.replace(".bam", ".unique.bam")
    sp.call(rmdup_cmd, shell = True)
    sp.call("samtools index " + filename.replace(".bam", ".unique.bam"), shell = True)



def dedup_bams(bams_dir, n_cores):
    bam_files = os.listdir(bams_dir)
    bam_files = [bams_dir + "/" + x for x in bam_files]

    pool = mp.Pool(n_cores)
    pool.map(_dedup_index_bam_, bam_files)
    pool.close()



def unmerge_bams(out_file_prefix, valid_cells, n_cores):
    input_file = out_file_prefix + ".merged.aligned.sorted.bam"
    
    if not os.path.isdir(out_file_prefix + "_deindexed_bam"):
        os.mkdir(out_file_prefix + "_deindexed_bam")

    with ps.AlignmentFile(input_file, "rb") as f_in:
        indexed = ps.IndexedReads(f_in)
        indexed.build()
        header = f_in.header.copy()

        for cell_id in valid_cells:
            out_file_name = out_file_prefix + "_deindexed_bam/" + cell_id + ".bam"
            with ps.AlignmentFile(out_file_name, "wb", header = header) as f_out:
                iterator = indexed.find(cell_id)
                for i in iterator:
                    f_out.write(i)



def run_alignment(out_file_prefix, ref_path, n_cores):
    align_cmd = "bwa mem -t " + str(n_cores) + " " + ref_path + " " + \
        out_file_prefix + ".merged.fastq > " + out_file_prefix + ".merged.aligned.sam"
    sp.call(align_cmd, shell = True)
    
    sort_cmd = "samtools sort -@ " + str(n_cores) + " -o " + out_file_prefix + ".merged.aligned.sorted.bam " + \
        out_file_prefix + ".merged.aligned.sam"
    sp.call(sort_cmd, shell = True)
    sp.call("samtools index " + out_file_prefix + ".merged.aligned.sorted.bam", shell = True)
    sp.call("rm " + out_file_prefix + ".merged.aligned.sam", shell = True)


	
def merge_fastq(fastq_dir, out_file_prefix):
    read_counts = Counter()
    out_file = out_file_prefix + ".merged.fastq"
    with open(out_file, "w") as f_out:
        for fi in os.listdir(fastq_dir):
            fi_name = fi.replace("_R1.fastq", "")
            with ps.FastxFile(fastq_dir + "/" + fi) as f_in:
                for read in f_in:
                    read_counts[fi_name] += 1
                    read.name = fi_name
                    f_out.write(str(read) + "\n")

    return read_counts
    


def run_deindexer(out_file_prefix, read1_fastq_file, i7_index_list, i5_index_list, r5_index_list,
		  n_mismatch = 1):
    if not os.path.isdir(out_file_prefix + "_deindexed_fastq"):
        os.mkdir(out_file_prefix + "_deindexed_fastq")
    
    deindexer_cmd = "deindexer -f rbbb -c scTHS.cfg.txt -b " + i7_index_list + " -b " + \
        i5_index_list + " -b " + r5_index_list + " " + read1_fastq_file + " " + \
        out_file_prefix + ".i7.idx.fastq " + out_file_prefix + ".i5.idx.fastq " + \
        out_file_prefix + ".r5.idx.fastq -o " + out_file_prefix + "_deindexed_fastq" + \
        " -m " + str(n_mismatch)

    sp.call("ulimit -n 500000; " + deindexer_cmd, shell = True)


	
def build_config(i7_cfg_file, i5_cfg_file, r5_cfg_file, out_file_prefix):
    i7_indexes = []
    with open(i7_cfg_file) as f:
        for line in f:
            i7_indexes.append(str.strip(line))
    
    i5_indexes = []
    with open(i5_cfg_file) as f:
        for line in f:
            i5_indexes.append(str.strip(line))

    r5_indexes = []
    with open(r5_cfg_file) as f:
        for line in f:
            r5_indexes.append(str.strip(line))

    indexes_list = []
    indexes_list.append(list(range(0,len(i7_indexes))))
    indexes_list.append(list(range(0,len(i5_indexes))))
    indexes_list.append(list(range(0,len(r5_indexes))))
    
    combined_indexes = list(product(*indexes_list))
    
    with open(out_file_prefix + ".cfg.txt", "w") as f:
        for prod in combined_indexes:
            i,j,k = prod
            out_line = str(i + 1) + "\t" + str(j + 1) + "\t" + str(k + 1) + "\t" + \
                    i7_indexes[i] + "_" + i5_indexes[j] + "_" + r5_indexes[k] + "\n"
            f.write(out_line)
    

	
def split_index(idx_file, out_file_prefix, start_idx, stop_idx):
    f_in = ps.FastxFile(idx_file)
    f_out = open(out_file_prefix + ".idx.fastq", mode = "w")
    for read in f_in:
        read.sequence = read.sequence[start_idx:stop_idx]
        read.quality = read.quality[start_idx:stop_idx]
        f_out.write(str(read) + "\n")
    f_in.close()
    f_out.close()



if __name__ == "__main__": main()
