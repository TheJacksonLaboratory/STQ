# Clear err_jupyter.err before starting if it exists
[ -f err_jupyter.err ] && > err_jupyter.err
sbatch run-jupyter-notebook.sb

for i in $(seq 1 60); do
    printf "\r%2ds " "$i"
    sleep 1
    url=$(cat err_jupyter.err 2>/dev/null | grep "http" | grep -v "]" | grep -v "127" | tail -1)
    if [ -n "$url" ]; then
        printf "\r%s\n" "$url"
        break
    fi
done
[ -z "$url" ] && printf "\rServer could not start\n"
