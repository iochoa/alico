This README file will detail the usage and purpose of the SAM compressor AliCo.

AliCo is a compressor and decompressor for SAM files and reference files. It can be downloaded from github.com/iochoa/Streamer_SAMfiles.

To make the program, download the code from GitHub. Make sure to have a version of cmake 3.5.0 or more recent installed and loaded. Make sure to be in the directory "Streamer_SAMfiles". The best practice for cmake is to create an empty directory "build-dir", change to that directory, and run "cmake .." to create the Makefile in "Streamer_SAMfiles". When that is done, change to the parent directory and simply run "make" to create the executable program SAMStreamer.

Running "./SAMStreamer -h" will bring up this help message:

Usage: ./SAMStreamer (options) [input file] [output] [ref file]
Options are:
	-x		: Regenerate file from compressed file
	-c [ratio]	: Compress using [ratio] bits per bit of input entropy per symbol
	-d [rate]	: Decompress and Download file from remote [output]
	-u [rate]	: Compress and Upload file to remote [output]
	-s [rate]	: Stream file to remote [output]
	-D [M|L|A]	: Optimize for MSE, Log(1+L1), L1 distortions, respectively (default: MSE)
	-q		: CALQ Mode
	-p [num]	: Set the polyploidy value for CALQ (default: 2)
	-a [num]	: Set the max quality value for CALQ (default: 41)
	-i [num]	: Set the min quality value for CALQ (default: 0)
	-o [num]	: Set the quality value offset for CALQ (default: 33)
	-e		: Compress reference file along with SAM file
	-m		: Decompress reference file along with SAM file
	-h		: Print this help

The default is compression mode (-c), so this need not be specified for lossless compression.
 
For normal compression mode, run "./SAMStreamer <SAM file name> <output directory> <reference file name>". This will create a few files in the output directory. "headers" contains the compressed headers and "mapped_reads" contains the compressed mapped reads from the SAM file. "unmapped_reads.gz" is a file uncompressed by the SAMStreamer but compressed by gzip.

For normal decompression mode, run "./SAMStreamer <input directory> <regenerated SAM file name> <reference file name>". This reads the directory automatically for the correct files and writes the final decompressed file to the specified location.

OPTIONS
    -c [ratio], compression
        Specify the bits per bit of input entropy per symbol in compression mode. 

    -u, upload
        Compress and upload the file to remote.

    -s, stream
        Stream the file to remote.

    -D [type], distortion optimization
        Specify the distortion optimization.  The types are:
            M : MSE distortions
            L : Log(1+L1) distortions
            A : L1 distortions
            If not specified, the default is MSE.

    -x, decompress
        Decompress a compressed SAM file.

    -e, compress reference file
        Compress the reference file as well as the SAM file. This creates the additional files "reference_num" and "reference_comp" which contain compressed information about the reference file.

    -m, decompress reference file
        Decompress the reference file which has been compressed by the above method. There is no reference file argument required in this mode, simply run "./SAMStreamer -m <input directory> <regenerated SAM file name>". This will create the regenerated reference file under the name "reference_local.fa" in the input directory specified.

CALQ OPTIONS

    -q, CALQ compression
        Compress the quality values with CALQ. This creates an additional file in the output directory "quality_values_calq" which contains the compressed quality values.

    -z, CALQ Decompression
        Decompress the quality values with CALQ. Usage: "./SAMStreamer -z <input Directory> <output file name> <reference file name>".

    -p, polyploidy
        Set the polyploidy value for CALQ (2 by default).

    -a, max quality value
        Set the maximum quality value for CALQ (41 by default).

    -i, min quality value
        Set the minimum quality value for CALQ (0 by default).

    -o, quality value offset
        Set the quality value offset for CALQ (33 by default).



