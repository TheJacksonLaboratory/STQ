version=(`nextflow -v`)
major=$(echo ${version[2]} | cut -d. -f1)

if [ "$major" -lt "24" ]; then
    read -p "Update nextflow to use the pipeline. Proceed? (y/n): " confirm && [[ $confirm == [yY] ]] && nextflow self-update
fi