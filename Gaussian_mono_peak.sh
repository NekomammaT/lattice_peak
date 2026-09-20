for ((i=0;i<5;i++))
do
    OMP_NUM_THREADS=4 ./Gaussian_mono_peak $i
done
