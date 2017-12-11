#/bin/bash
if [[ -z "$1" ]]; then
    echo "No action supplied"
elif [[ "$1" == "to" ]]; then
  echo "Upload source code to" $2
  if [[ "$2" == "lanl" ]]; then
    scp -r Makefile src qsub apps lanl:ar-tn.lanl.gov:/usr/projects/cint/cint_sces/tensorlib/apps/
  else
    rsync -avzh \
    --include-from 'code-include.txt' \
    --exclude-from 'code-exclude.txt' \
    $(pwd) $2:~/GitRepo/
  fi
elif [[ "$1" == "from" ]]; then
  echo "Get data from" $2
  if [[ "$2" == "lanl" ]]; then
    scp -r lanl:ar-tn.lanl.gov:/net/scratch3/chengyanlai/$3 data/
  else
    rsync -avzh --update \
    --exclude-from 'data-exclude.txt' \
    $2:~/GitRepo/ScaLapackSVD/data/$3 data/
  fi
fi
