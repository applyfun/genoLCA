biaslevels=(bias_level_1 bias_level_2 bias_level_3)

for i in "${biaslevels[@]}"; do
	cp simulation-template.sh j-${i}.sh
	sed -i "s/IDX/${i}/g" j-${i}.sh
done


