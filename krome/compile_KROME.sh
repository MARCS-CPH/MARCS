network=networks/non_eq/react_Chapman_photo_krome
project=react_Chapman_photo_krome
options=options_photochem

if [ -d build_$project ] ; then #delete old build folder
rm build_$project/*
fi 

./krome -n $network -project=$project -options=$options

if [ -d MARCS_build ] ; then #delete old MARCS_build
echo "Deleting old MARCS_build folder"
rm -r MARCS_build
fi

if [ ! -d MARCS_build ] ; then #create new MARCS_build folder if none exists
echo "Creating new MARCS_build folder"
mkdir MARCS_build
fi

cp build_$project/* MARCS_build #copy over files from the krome build folder
mv $project.kpj MARCS_build
echo "Copy over files to MARCS_build"
