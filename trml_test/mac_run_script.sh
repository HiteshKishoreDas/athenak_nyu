export athenak=/Users/hitesh/hitesh/git/athenak_nyu
export build=$athenak/build

export OMP_NUM_THREADS=4

$build/src/athena -i $athenak/trml_test/athinput.TRML_simple -d $athenak/trml_test
