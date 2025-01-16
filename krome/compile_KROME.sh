project=react_Chapman_incl_photo
network=networks/non_eq/react_Chapman_incl_photo
options=options_standard

if [ -d build_$project ] ; then
rm build_$project/*
fi 

./krome -n $network -project=$project -options=$options

if [ ! -d MARCS_build ] ; then
mkdir MARCS_build
fi

cp build_$project/* MARCS_build
mv $project.kpj MARCS_build
cp MARCS_build/reactions_verbatim.dat ../
