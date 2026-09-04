#!/usr/bin/fish

cd cloudy_grid;
for run in run_*
    echo $run;
    cd $run/skirt;
    skirt model.ski;
    cd ../..;
end