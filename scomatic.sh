cd /data/AMLscRNA/1503/outs
output_dir=/data/AMLscRNA/1503/outs/scomatic/step1
SCOMATIC=/data/hello_miniconda/SComatic
head -1 /data/AMLscRNA/1026/outs/p1b.txt>/data/AMLscRNA/1503/outs/P6a.txt
cat /data/AMLscRNA/mergeblood/20220501/meta.bed|cut -f1,26|grep 'P6t'|sed 's/P6t_//g'|sed 's/-1//g'>>/data/AMLscRNA/1503/outs/P6a.txt
/data/hello_miniconda/home/envs/SComatic/bin/python /data/hello_miniconda/SComatic/scripts/SplitBam/SplitBamCellTypes.py --bam possorted_genome_bam.bam \
        --meta P6a.txt \
        --id P6a \
        --n_trim 5 \
        --max_nM 5 \
        --max_NH 1 \
        --outdir $output_dir

REF=/data/ref/hg38.fa
output_dir1=/data/AMLscRNA/1503/outs/scomatic/step1
output_dir2=/data/AMLscRNA/1503/outs/scomatic/step2

for bam in $(ls -d $output_dir1/*bam);do
  
  # Cell type
  cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')

  # Temp folder
  temp=$output_dir2/temp_${cell_type}
  mkdir -p $temp

  # Command line to submit to cluster
  /data/hello_miniconda/home/envs/SComatic/bin/python /data/hello_miniconda/SComatic/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
    --ref $REF \
    --chrom all \
    --out_folder $output_dir2 \
    --min_bq 30 \
    --tmp_dir $temp \
    --nprocs 20

  rm -rf $temp
done

sample=P6a
output_dir3=/data/AMLscRNA/1503/outs/scomatic/step3
output_dir2=/data/AMLscRNA/1503/outs/scomatic/step2
/data/hello_miniconda/home/envs/SComatic/bin/python /data/hello_miniconda/SComatic/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder ${output_dir2} \
  --outfile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv

output_dir4=/data/AMLscRNA/1503/outs/scomatic/step4

sample=P6a
REF=/data/ref/hg38.fa

/data/hello_miniconda/home/envs/SComatic/bin/python /data/hello_miniconda/SComatic/scripts/BaseCellCalling/BaseCellCalling.step1.py \
          --infile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv \
          --outfile ${output_dir4}/${sample} \
          --ref $REF

# Step 4.2
editing=/data/hello_miniconda/SComatic/RNAediting/AllEditingSites.hg38.txt
PON=/data/hello_miniconda/SComatic/PoNs/PoN.scRNAseq.hg38.tsv

/data/hello_miniconda/home/envs/SComatic/bin/python /data/hello_miniconda/SComatic/scripts/BaseCellCalling/BaseCellCalling.step2.py \
          --infile ${output_dir4}/${sample}.calling.step1.tsv \
          --outfile ${output_dir4}/${sample} \
          --editing $editing \
          --pon $PON

bedtools intersect -header -a ${output_dir4}/${sample}.calling.step2.tsv -b /data/hello_miniconda/SComatic/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed | awk '$1 ~ /^#/ || $6 == "PASS"' > ${output_dir4}/${sample}.calling.step2.pass.tsv
