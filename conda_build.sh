if [ -z ${CONDA_BUILD+x} ]; then
	source /opt/conda/conda-bld/dfast_1563057817011/work/build_env_setup.sh
fi
#!/bin/sh

APPROOT=$PREFIX/opt/$PKG_NAME-$PKG_VERSION

rm bin/*/aragorn
rm bin/*/barrnap
rm -r bin/Darwin/barrnap-0.8
rm bin/*/blastdbcmd
rm bin/*/blastp
rm bin/*/ghost*
rm bin/*/hmm*
rm -r bin/Darwin/lib
rm bin/*/makeblastdb
rm bin/*/mga
rm bin/*/rpsblast


mkdir -p $APPROOT
cp -r ./* $APPROOT

ln -s ${APPROOT}/dfast ${PREFIX}/bin/
ln -s ${APPROOT}/scripts/file_downloader.py ${PREFIX}/bin/dfast_file_downloader.py
ln -s ${APPROOT}/scripts/file_downloader.py ${PREFIX}/bin
