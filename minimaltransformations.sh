# brain extraction
parallel --dry-run --plus --bar fslmaths {} -bin {..}_mask ::: /path/to/*brain_1mm.nii.gz
# B1 field bias correction
parallel --dry-run --plus N4BiasFieldCorrection -d 3 -i {} -x {..}_mask.nii.gz -r -o {..}_biasfieldcorrected.nii.gz ::: /path/to/*brain_1mm.nii.gz
# Linear registration to MNI
ls /path/to/*brain_1mm_biasfieldcorrected.nii.gz | parallel --dry-run --plus flirt -in {} -ref $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz -cost normmi -omat {..}_flirt.mat -out {..}_intemplate.nii.gz