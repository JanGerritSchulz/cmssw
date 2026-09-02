echo '#### Test Writing and Reading ZVertexSoA'

scriptdir=$CMSSW_BASE/src/DataFormats/VertexSoA/test/

echo '> Writing'

cmsRun ${scriptdir}/testWriteHostZVertexSoA.py testZVertexSoA.root

if [ $? -ne 0 ]; then
    exit 1;
fi

echo '> Reading'

cmsRun ${scriptdir}/testReadHostZVertexSoA.py testZVertexSoA.root

if [ $? -ne 0 ]; then
    exit 1;
fi

echo '>>>> Done! <<<<'

