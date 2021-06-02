# Examples, tutorials

Currently the file [210427.campy_bordetella.ipynb](210427.campy_bordetella.ipynb) contains the `jupyter` notebook for the
*Campylobacter* and *Bordetella* homopolymer tract analyses conducted with tatajuba, showing how to summarise and
display the output information. Since the output and extra files are too big, they are compacted using `xz` into file
`210427.outdir.txz`. To decompress it, run
```
tar Jxvf 210427.outdir.txz
```
It will create two directories and consume 400MB of disk space. The notebook expects `gzip`ed files, you may need to compress the `tsv` files or change the notebook. 
