project=react_Chapman_incl_photo
if [ -d build_$project ] ; then
rm build_$project/*
fi 

./krome -n networks/non_eq/react_Chapman_incl_photo -useN -noSinkCheck -checkConserv -project=$project -noExample -unsafe -ATOL=1d-40 -RTOL=1D-8

if [ ! -d MARCS_build ] ; then
mkdir MARCS_build
fi

cp build_$project/* MARCS_build
mv $project.kpj MARCS_build
cp MARCS_build/reactions_verbatim.dat ../
