project=MARCS_NO_transport
#./krome -n networks/react_NO -useN -reverse -noSinkCheck -checkConserv -project=$project -noExample -unsafe -ATOL=1d-30 -RTOL=1D-6
./krome -n networks/react_NO -useN -reverse -noSinkCheck -project=$project -noExample -unsafe -ATOL=1d-30 -RTOL=1D-6
#./krome -n networks/react_NO -useN -reverse -noSinkCheck -project=$project -noExample -unsafe
#./krome -n networks/react_NO -useN -reverse -noSinkCheck -checkConserv -project=$project -noExample -unsafe -ATOL=1d-30 -RTOL=1D-6 -useEquilibrium
#./krome -n networks/react_NO_transport -useN -noSinkCheck -project=$project -noExample -unsafe -nomassCheck -ATOL=1d-30 -RTOL=1D-6

if [ ! -d MARCS_build ] ; then
mkdir MARCS_build
fi

cp build_$project/* MARCS_build
mv $project.kpj MARCS_build
cp MARCS_build/reactions_verbatim.dat ../
