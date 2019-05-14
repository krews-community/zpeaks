package step

import htsjdk.samtools.*
import model.*
import mu.KotlinLogging
import java.nio.file.Path


private val log = KotlinLogging.logger {}

enum class PileUpAlgorithm {
    // Use only the start of each alignment
    START,
    // Use the mid-point between start and end for each alignment. *Paired End Only
    MID_POINT,
    // Use every BP between start and end for each alignment. *Paired End Only
    LENGTH
}

// When using the LENGTH pile-up algorithm we will not include values for lengths that are larger than this value.
// They are most likely junk and will waste resources.
const val LENGTH_LIMIT = 100_000

class PileUp(
    // Pile-up values in chromosome.
    private val values: IntArray,
    // chromosome length pulled directly from BAM File
    val chromosomeLength: Int,
    // Sum of all pile-up values in chromosome calculated on the fly and cached for efficiency
    val sum: Int
) {
    operator fun get(bp: Int): Int = values[bp]
}

/**
 * Reads a SAM or BAM file and creates a pile-up representation in memory.
 *
 * @param samPath: Path to the SAM or BAM file
 * @param strand: Strand that we want to count using pile-up algorithm
 * @param atacMode: Whether or not this alignment is for an ATAC experiment. Algorithm will use special offsets
 * if true.
 * @param pileUpAlgorithm: Algorithm we use to choose the values that we pile up.
 * @return the in-memory pile-up
 */
fun pileUpSam(samPath: Path, strand: Strand, atacMode: Boolean, pileUpAlgorithm: PileUpAlgorithm): Map<String, PileUp> {
    val samReader = SamReaderFactory.make().open(samPath)
    val chromosomeLengths: Map<String, Int> = samReader.fileHeader.sequenceDictionary.sequences
        .map { it.sequenceName to it.sequenceLength }.toMap()
    val values = mutableMapOf<String, IntArray>()
    val sums = mutableMapOf<String, Int>()

    // use block to auto-close
    samReader.use { reader ->
        reader.forEach { record ->
            // If we're only using plus strand values and this record is a plus strand and
            // visa versa for minus strand, continue.
            if ((strand == Strand.PLUS && record.readNegativeStrandFlag) ||
                (strand == Strand.MINUS && !record.readNegativeStrandFlag)) return@forEach

            // Chromosome name / key
            val chr = record.referenceName

            if (!values.containsKey(chr)) values[chr] = IntArray(chromosomeLengths.getValue(chr)) { 0 }
            sums.putIfAbsent(chr, 0)

            val chrValues = values.getValue(chr)
            val chrLength = chromosomeLengths.getValue(chr)
            val start = pileUpStart(record, atacMode, chrLength)
            when (pileUpAlgorithm) {
                PileUpAlgorithm.START -> {
                    chrValues[start]++
                    sums[chr] = sums.getValue(chr) + 1
                }
                PileUpAlgorithm.MID_POINT -> {
                    val end = pileUpEnd(record, atacMode, chrLength)
                    val midPoint = (start + end) / 2
                    chrValues[midPoint]++
                    sums[chr] = sums.getValue(chr) + 1
                }
                PileUpAlgorithm.LENGTH -> {
                    val end = pileUpEnd(record, atacMode, chrLength)
                    val length = end - start
                    if (length > LENGTH_LIMIT) return@forEach
                    for (i in start until end) {
                        chrValues[i]++
                    }
                    sums[chr] = sums.getValue(chr) + length
                }
            }
        }
    }

    return values.keys.map { chr ->
        chr to PileUp(
            chromosomeLength = chromosomeLengths.getValue(chr),
            values = values.getValue(chr),
            sum = sums.getValue(chr)
        )
    }.toMap()
}

/**
 * Calculate the "start" value for our pile-up algorithm for the given record and settings.
 *
 * This start will actually be the record's end it's strand is minus
 */
private fun pileUpStart(record: SAMRecord, atacMode: Boolean, chrLength: Int): Int {
    return if (!record.readNegativeStrandFlag) {
        if (atacMode) atacStart(record, chrLength) else record.start
    } else {
        if (atacMode) atacEnd(record) else record.end
    }
}

/**
 * Calculate the "end" value for our pile-up algorithm for the given SAM record and settings.
 *
 * This end will actually be the record's start it's strand is minus
 */
private fun pileUpEnd(record: SAMRecord, atacMode: Boolean, chrLength: Int): Int {
    return if (!record.readNegativeStrandFlag) {
        if (atacMode) atacEnd(record) else record.end
    } else {
        if (atacMode) atacStart(record, chrLength) else record.start
    }
}

/**
 * Calculate a SAM record's start adjusted for ATAC Mode
 */
private fun atacStart(record: SAMRecord, chrLength: Int): Int =
    if (record.start + 4 > chrLength) chrLength else record.start + 4

/**
 * Calculate a SAM record's end adjusted for ATAC Mode
 */
private fun atacEnd(record: SAMRecord): Int = if (record.end - 5 < 0) 0 else record.end - 5
