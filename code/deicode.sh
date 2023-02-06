qun_dir=/mnt/zarrinpar/scratch/sfloresr/HE/HE_metab/liverpath/data_files/LCMS/quantification_table/
mtb_dir=/mnt/zarrinpar/scratch/sfloresr/HE/HE_metab/liverpath/data_files/LCMS/DB_result/
md_dir=/mnt/zarrinpar/scratch/sfloresr/HE/HE_metab/liverpath/data_files/LCMS/metadata_table/
outpt_dir=/mnt/zarrinpar/scratch/sfloresr/HE/HE_metab/liverpath/diversity_metrics_LCMS/

#deicode calculation of TIPS progression for the hepatic vein samples

biom convert \
  -i $qun_dir/quantification_table-hepatic.txt \
  -o $qun_dir/quantification_table-hepatic.biom \
  --table-type="OTU table" \
  --to-hdf5

qiime tools import \
  --input-path $qun_dir/quantification_table-hepatic.biom \
  --output-path $qun_dir/quantification_table-hepatic.qza \
  --type FeatureTable[Frequency]

qiime deicode rpca \
    --i-table $qun_dir/quantification_table-hepatic.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot $outpt_dir/ordination-hepatic.qza \
    --o-distance-matrix $outpt_dir/beta_deicode_lcms_hepatic.qza

qiime diversity beta-group-significance \
    --i-distance-matrix $outpt_dir/beta_deicode_lcms_hepatic.qza \
    --m-metadata-file $md_dir/metadata_table-hepatic.txt \
    --m-metadata-column ATTRIBUTE_blood_procedure \
    --p-method permanova \
    --p-pairwise \
    --o-visualization $outpt_dir/beta_deicode_lcms_hepatic-significance.qzv


#deicode calculation of TIPS progression for the peripheral vein samples

biom convert \
  -i $qun_dir/quantification_table-peripheral.txt \
  -o $qun_dir/quantification_table-peripheral.biom \
  --table-type="OTU table" \
  --to-hdf5

qiime tools import \
  --input-path $qun_dir/quantification_table-peripheral.biom \
  --output-path $qun_dir/quantification_table-peripheral.qza \
  --type FeatureTable[Frequency]

qiime deicode rpca \
    --i-table $qun_dir/quantification_table-peripheral.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot $outpt_dir/ordination-peripheral.qza \
    --o-distance-matrix $outpt_dir/beta_deicode_lcms_peripheral.qza

qiime diversity beta-group-significance \
    --i-distance-matrix $outpt_dir/beta_deicode_lcms_peripheral.qza \
    --m-metadata-file $md_dir/metadata_table-peripheral.txt \
    --m-metadata-column ATTRIBUTE_blood_procedure \
    --p-method permanova \
    --p-pairwise \
    --o-visualization $outpt_dir/beta_deicode_lcms_peripheral-significance.qzv

#deicode calculation of TIPS progression for the peripheral vein samples (just pre vs. post)

biom convert \
  -i $qun_dir/quantification_table-peripheral-prepost.txt \
  -o $qun_dir/quantification_table-peripheral-prepost.biom \
  --table-type="OTU table" \
  --to-hdf5

qiime tools import \
  --input-path $qun_dir/quantification_table-peripheral-prepost.biom \
  --output-path $qun_dir/quantification_table-peripheral-prepost.qza \
  --type FeatureTable[Frequency]

qiime deicode rpca \
    --i-table $qun_dir/quantification_table-peripheral-prepost.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot $outpt_dir/ordination-peripheral-prepost.qza \
    --o-distance-matrix $outpt_dir/beta_deicode_lcms_peripheral-prepost.qza

#deicode calculation of HE grade for the post hepatic vein samples

biom convert \
  -i $qun_dir/quantification_table-post-hepatic.txt \
  -o $qun_dir/quantification_table-post-hepatic.biom \
  --table-type="OTU table" \
  --to-hdf5

qiime tools import \
  --input-path $qun_dir/quantification_table-post-hepatic.biom \
  --output-path $qun_dir/quantification_table-post-hepatic.qza \
  --type FeatureTable[Frequency]

qiime deicode rpca \
    --i-table $qun_dir/quantification_table-post-hepatic.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot $outpt_dir/ordination-post-hepatic.qza \
    --o-distance-matrix $outpt_dir/beta_deicode_lcms_post-hepatic.qza

qiime diversity beta-group-significance \
    --i-distance-matrix $outpt_dir/beta_deicode_lcms_post-hepatic.qza \
    --m-metadata-file $md_dir/metadata_table-post-hepatic.txt \
    --m-metadata-column WPostTIPS_HE \
    --p-method permanova \
    --p-pairwise \
    --o-visualization $outpt_dir/beta_deicode_lcms_post-significance-hepatic.qzv

#run qurro on the deicode results of TIPS progression in the peripheral vein samples

qiime qurro loading-plot \
    --i-table $qun_dir/quantification_table-peripheral.qza \
    --i-ranks $outpt_dir/ordination-peripheral.qza \
    --m-sample-metadata-file $md_dir/metadata_table-peripheral.txt \
    --m-feature-metadata-file $mtb_dir/DB_result/spectral_matches_clean_v2.txt \
    --o-visualization $outpt_dir/qurro_deicode_peripheral_v2.qzv
