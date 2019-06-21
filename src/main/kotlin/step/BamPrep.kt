package step

import htsjdk.samtools.*
import java.nio.file.Path

/**
 * Check each bam for the existence of an index. If it does not exist, it will be created.
 *
 * @param bams: bam files.
 * @param chrFilter: List of chromosomes we want to process and therefor will return the lengths for.
 *
 * @returns a list of chromosomes with max lengths (The highest length given in a bam per chromosome)
 */
fun prepBams(bams: List<Path>, chrFilter: List<String>?): Map<String, Int> {
    val chrLengths = mutableMapOf<String, Int>()
    for(bam in bams) {
        val samReader = SamReaderFactory.makeDefault()
            .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
            .open(bam)
        samReader.use { reader ->
            if (SamFiles.findIndex(bam) == null) {
                BAMIndexer.createIndex(reader, bam.resolveSibling("${bam.fileName}.bai"))
            }
            for(sequence in reader.fileHeader.sequenceDictionary.sequences) {
                val chrName = sequence.sequenceName
                val chrLength = sequence.sequenceLength
                if (chrFilter != null && !chrFilter.contains(chrName)) continue
                if(!chrLengths.containsKey(chrName) || chrLengths.getValue(chrName) > chrLength) {
                    chrLengths[chrName] = chrLength
                }
            }
        }
    }
    return chrLengths
}