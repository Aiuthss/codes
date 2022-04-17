# $ID_DATA="GSE194380"
# $STRAIN="MMusculus"
$ID_DATA = Read-Host "IDを入力してください"
$STRAIN = Read-Host "Strainを入力してください"

$NAME_IMAGE="rna-seq"
$NAME_CONTAINER="RNA-seq_automation"

#docker build -t $NAME_IMAGE C:\RNA-seq_data\codes

docker run -v C:\RNA-seq_data\:/home/C -v D:\RNA-seq_archives:/home/D --name $NAME_CONTAINER -dit $NAME_IMAGE

#ID_CONTAINER=`docker container ls -q -f name=$NAME_CONTAINER`

docker exec -d ${NAME_CONTAINER} /bin/bash -c "bash /home/C/codes/run.sh ${ID_DATA} ${STRAIN} && Rscript /home/C/codes/summarize.R ${ID_DATA} ${STRAIN}"
#docker stop RNA-seq_automation
#docker container rm RNA-seq_automation