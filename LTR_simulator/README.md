# LTR_simulator
A tool that simulates LTR insertions (intact and/or solo) and mutations. 

`nohup bash ./LTR_simulator/LTR_simulator.sh -sL ./all.pass.fa -sLc 3000 -h 5000000 -l 500000 -s 500000 -uL ./all.young.pass.fa &`  

# Extracting LTRs for use as input.
`cat *_mod.pass.list | grep 'pass' | awk '$8 ~ /1.0000/' | sort | uniq > all.young.pass.list `  
`cat *_mod.pass.list | grep 'pass' | sort | uniq > all.pass.list `  
`cat all.young.pass.list | cut -f 1 | awk 'BEGIN {FS=OFS="\t"} {if (match($1, /(.+):([0-9]+)..([0-9]+)/, arr)) {id=arr[1]; start=arr[2]; end=arr[3]; if (start > end) {tmp=start; start=end; end=tmp} $1=id ":" start ".." end} print}' | sed -e 's/:/\t/g' -e 's/\.\./\t/g' > all.young.pass.bed `  
`cat all.pass.list | cut -f 1 | awk 'BEGIN {FS=OFS="\t"} {if (match($1, /(.+):([0-9]+)..([0-9]+)/, arr)) {id=arr[1]; start=arr[2]; end=arr[3]; if (start > end) {tmp=start; start=end; end=tmp} $1=id ":" start ".." end} print}' | sed -e 's/:/\t/g' -e 's/\.\./\t/g' > all.pass.bed `  
`bedtools getfasta -fi all.fa -fo all.young.pass.fa -bed all.young.pass.bed `  
`bedtools getfasta -fi all.fa -fo all.pass.fa -bed all.pass.bed `  

### Example run.  
`bash LTR_simulator/LTR_simulator.sh -g1 9311.fa -g2 Nipponbare.fa -t 100 -lr 1e-12 -r 0:1 -l 1000000 -s 2000000 -h 10000000 ` 

- Fragmented TEs in a real genome add to off target TE annotations. For my purpose, its better to use a synthetic genome backbone to insert TEs into. This allows more predictable annotations and easier quantification. `synthetic_genome.py` can build a fake genome.
  - I may need to insert genes to guide syntenic anchoring. Currently, without genes, syentic anchors only cover a portion.

- Similarly, not all LTRs are predicted with equal reliability by LTR_retriever. I used the MTEC & rice curated libraries rather than LTR annotations generated in-house. I inserted those into a fake genome, annnotated with LTR_retriever, and retained only those curated LTRs that were recapitulated by LTR_retriever.

- `python LTR_simulator/bin/ltr_inserter.py -genome noTEs.fa -LTR LTR_simulator/curated_LTRs_proof.fa -n 800`
  - Insert 800 LTRs into the fake genome at random locations.
  - These are LTRs of any age.
  
- Flanking LTRs have a tag `CCCCCCCCCCCCCCCCCC`.
  - `cat pairwise_synLTR/N8ha1m.N8hb1m_mod.pass.list.count.no_dups.xls | grep 'Site_type:Fl' | cut -f 1 | sed -e 's/:/\t/g' -e 's/\.\./\t/g' > flank.bed`
  - `bedtools getfasta -fi pre_pairwise_evaluation/N8ha1m.N8hb1m.fa -fo flank.fa -bed flank.bed`
  - `cat flank.fa | grep 'CCCCCCCCCCCCCCCCCC' | wc -l`

### Steps.  
(1) Mutate the genome using a fixed rate, X, based on the mutation rate and number of generations.  
(2) Mutate the LTRs using a variable rate where X is max.  
(3) Insert the LTRs into the genome.  

### Outputs. 'l', 's', and 'h' control the number of generations simulated.  
`93111m.fa` # One million generations simulated.  
`93113m.fa`  # Three million generations simulated.  
`93115m.fa` # Five million generations simulated.  
`93117m.fa` # Seven million generations simulated.  
`93119m.fa` # Nine million generations simulated.  
`Nipponbare1m.fa` # One million generations simulated.  
`Nipponbare3m.fa` # Three million generations simulated.  
`Nipponbare5m.fa` # Five million generations simulated.  
`Nipponbare7m.fa` # Seven million generations simulated.  
`Nipponbare9m.fa` # Nine million generations simulated.  

### The tool doesnt simulate LTR removals.   
- For this, i'd need to either track where I had inserted the fake LTRs (may be tough due to coordinate shifting at each new insertion) or run LTR_retriever.   
  - If I had coordinates for all LTRs, I could write a new script to remove them.  
  - The LTR insertion rate of `1e-11` was selected as default.   
  - In reality, the rate is likely higher, but I found higher rates tend to give an unreasonable number of insertions.  
  - This is becasue we dont simulate LTR removals, so our insertion rate is actually the rate of retained insertions, which would be lower than the number of total insertions.   

### Ideas for improvement.  
(1) ltr_mutator with `mode 0` should print the number of mutations that it adds to the genome. This will allow me to more easily calculate rates.   
(2) Simulate LTR removals.   
(3) Remove intermediate files.  
