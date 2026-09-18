for ((i=2;i<6;i++))
do
    OMP_NUM_THREADS=4 ./Gaussian_mono_peak $i
done
