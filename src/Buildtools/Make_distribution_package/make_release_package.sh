#this script creates a distribution package from the current working copy
#this script has to be called in the root of the working copy, i.e. in the directory containing src, Docs,...
#the first commandline argument is the destination path of the packed archive, e.g. ~/CPPAW_distrib.tar.gz
#the second commandline argument is the version identifier
#the third commandline argument is the path of the parms-file to be included, e.g. /home/bla/parms.blub_bla
# multiple parms-files are possible, the fist one is used to create docs, so it should be working on the current machine
#as temporary storage $WD is used 
#example: sh src/Buildtools/Make_distribution_package/make_distribution_package.sh /home/user0/tmp/cppaw_distrib.tar.gz v2023.1 /home/user0/Data/git/parms.linux

OUT="$1"
VERSION="$2"
PARMS="$3"
echo $OUT
WD="/tmp/CPPAW_distrib"

rm $OUT
rm -rf $WD
mkdir -p $WD

cp -r * $WD
cp -r .git $WD

cp -v $PARMS $WD

cd $WD
sh src/Buildtools/Version/getversion.sh
#replace DEVELOPMENT VERSION with 
dat=`tail -n 1 src/version.info`

echo "RELEASE VERSION" > src/version.info
echo "$VERSION" >>  src/version.info
echo "$dat" >> src/version.info
echo "https://github.com/cp-paw/cp-paw" >> src/version.info
#sh src/Buildtools/Version/getversion.sh src/version.info

rm -rf .git

./configure --with-parmfile=`basename $PARMS`
make docs
make clean
cd ..

tar -cvzf $OUT `basename $WD`
rm -rf $WD
