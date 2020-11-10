biaslevels=("bias_level_1" "bias_level_2" "bias_level_3")

for i in "${biaslevels[@]}"; do
	sbatch j-${i}.sh
	sleep 3
done