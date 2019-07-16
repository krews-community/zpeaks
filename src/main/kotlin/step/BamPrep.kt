package step

import htsjdk.samtools.*
import model.ChromBounds
import mu.KotlinLogging
import java.nio.file.Path

private val log = KotlinLogging.logger {}

/**
 * Check each bam for the existence of an index. If it does not exist, it will be created.
 *
 * @param bams: bam files.
 * @param chrFilter: List of chromosomes we want to process and therefor will return the lengths for.
 *
 * @returns a list of chromosomes with max lengths (The highest length given in a bam per chromosome)
 */
fun prepBams(bams: List<Path>, chrFilter: Map<String, IntRange?>? = null): Map<String, ChromBounds> {
    log.info { "Collecting chromosome metadata and checking BAMs for indexes..." }
    val chrLengths = mutableMapOf<String, Int>()
    for(bam in bams) {
        val samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
            .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
            .open(bam)
        samReader.use { reader ->
            if (SamFiles.findIndex(bam) == null) {
                log.info { "Index not found for ${bam.fileName}. Creating..." }
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

    val chrRanges = chrLengths.map { (chr, length) ->
        val range =
            if (chrFilter?.get(chr) != null) chrFilter.getValue(chr)!!
            else (0 until length)
        chr to ChromBounds(chr, length, range)
    }.toMap()

    log.info { "BAM prep complete!" }
    return chrRanges
}