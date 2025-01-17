echo $2
~/tools/timeout -m 4000000 ~/tools/dap/dap_src/dap-g \
    -d input/$1/$2.sbams.dat \
    -o output/$1/$2.dat \
    -l log/$1/$2.log \
    -msize 5 -ld_control 0.75 -size_limit 500
