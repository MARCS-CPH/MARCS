network=networks/non_eq/react_Chapman_photo_krome
project=react_Chapman_photo_krome
options=options_photochem

if [ -d build_$project ] ; then
rm build_$project/*
fi 

./krome -n $network -project=$project -options=$options

if [ ! -d MARCS_build ] ; then
mkdir MARCS_build
fi

cp build_$project/* MARCS_build
cp build_$project/* build
mv $project.kpj MARCS_build
cp MARCS_build/reactions_verbatim.dat ../
