# kim-ddSeq-barcode
BioRad_ddSeq_Barcodes. Project start: 9/11/17

To run script:

$ python parseBarcodes-N1.6.1.py -input ddSeq.sam -blocks barcodeBlocks.txt -output writeSamTest.sam

where:
parseBarcodes-N1.6.1.py is the driver script
ddSeq.sam is the SAM-formatted ddSeq sequence reads
barcodeBlocks.txt is the list of accepted Illumina barcode blocks
writeSamTest.sam is the output file with our barcode-containing reads
