for ((i=1;i<5;i++))
do
    OMP_NUM_THREADS=4 ./Gaussian_LN_peak $i
done
