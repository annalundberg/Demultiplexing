# Demultiplexing

This is an algorithm for demultiplexing sequence reads/indexes.

We need to determine if the indexes match per read. One sequence file is forward reads, the other is reverse. The reads appear in a consistent order between files. We need to filter the reads so that each index gets a file of forward reads and one of reverse reads. Since we don’t know which file has forward and which has reverse, we will arbitrarily assign. Reads with unclear indexes will be thrown into the pair of files with unknown indexes.


## Output: 
File names of all files created (2 files per index (forward and reverse) and 2 files of reads with unclear indexes), integer of properly matched indexes and float percentage of index hopping observed.'''
