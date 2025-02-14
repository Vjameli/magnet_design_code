#!/bin/csh

set DIR = ${HOME}/Robust/Thick_Matlab_5/test

foreach inst (`ls ${DIR}`)
    echo $DIR/$inst
    qsub -v "FICHIER=${DIR}/${inst}" script.csh

end
