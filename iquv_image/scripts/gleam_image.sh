#!/bin/bash

# dep1=6272402

# sleep 3600

for s in 1*
do
	# s=${f%.ms.zip}
	# mkdir $s
	# mv $f $s
	cat image.tmpl.sh | sed "s;OBSID;${s};g" > ${s}.img.sh
	# cat image_dep.tmpl.sh | sed "s;OBSID;${s};g" > ${s}.img_dep.sh
	mv ${s}.img.sh ${s}/
	# mv ${s}.img_dep.sh ${s}/
	cd ${s}/
	# XZ: I would like the successor job to start 10 min after the start of the precessor job.
	# tried --begin=now+10. The jobs were each delayed for 1 min!
	# dep=$(sbatch --depend=afterok:${dep1} ${s}.img.sh)
	# arrdep=(${dep// / })
	# dep1=${arrdep[3]}

	# dep=$(sbatch --depend=after:${dep1} ${s}.img_dep.sh)
	# arrdep=(${dep// / })
	# dep1=${arrdep[3]}

	sbatch ${s}.img.sh



	# echo $dep
	# echo $dep1
	cd ../
done
