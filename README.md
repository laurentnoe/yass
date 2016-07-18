# yass
(from  http://bioinfo.lifl.fr/yass/ ) 

YASS is a genomic similarity search tool, for nucleic (DNA/RNA) sequences in fasta or plain text format (it produces local pairwise alignments). Like most of the heuristic pairwise local alignment tools for DNA sequences (FASTA, BLAST, PATTERNHUNTER, BLASTZ/LASTZ, LAST ...), YASS uses seeds to detect potential similarity regions, and then tries to extend them to local alignments. This genomic search tool uses multiple transition constrained spaced seeds that enable to search more fuzzy repeats, as non-coding DNA/RNA. Another simple, but interesting feature is that you can specify the seed pattern used in the search step (as provided for example by [iedera](https://github.com/LaurentNoe/iedera))
