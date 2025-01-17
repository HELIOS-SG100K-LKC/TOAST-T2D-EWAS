ls -1 s1_bl_filter/ | parallel -j 20 --gnu --plus "cat s1_bl_filter/{} | awk '{if(\$5+\$6>=5) print \$0}' > s2_cov_gte5/{%.bismark.cov}_rmBL_cov5.bismark.cov;"
