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
    override val chrLength: Int,
    // Sum of all pile-up values in chromosome calculated on the fly and cached for efficiency
    val sum: Int
): SignalData {
    override operator fun get(bp: Int): Int = values[bp]
}

/**
 * @param strand: Strand that we want to count using pile-up algorithm
 * @param pileUpAlgorithm: Algorithm we use to choose the values that we pile up.
 * @param forwardShift: shifts forward (plus) strand by this amount
 * @param reverseShift: shifts reverse (minus) strand by this amount
 */
data class PileUpOptions(
    val strand: Strand,
    val pileUpAlgorithm: PileUpAlgorithm,
    val forwardShift: Int = 0,
    val reverseShift: Int = 0
)

/**
 * Reads a SAM or BAM file and creates a pile-up representation in memory.
 * @param samPath: Path to the SAM or BAM file
 * @return the in-memory pile-up
 */
fun runPileUp(samPath: Path, pileUpOptions: PileUpOptions): MutableMap<String, PileUp> {
    log.info { "Performing pile-up for sam file $samPath with strand=${pileUpOptions.strand}, " +
            "pileUpAlgorithm=${pileUpOptions.pileUpAlgorithm}, " +
            "forwardShift=${pileUpOptions.forwardShift}, reverseShift=${pileUpOptions.reverseShift}" }
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
            if ((pileUpOptions.strand == Strand.PLUS && record.readNegativeStrandFlag) ||
                (pileUpOptions.strand == Strand.MINUS && !record.readNegativeStrandFlag)) return@forEach

            // Chromosome name / key
            val chr = record.referenceName

            if (!values.containsKey(chr)) values[chr] = IntArray(chromosomeLengths.getValue(chr)) { 0 }
            sums.putIfAbsent(chr, 0)

            val chrValues = values.getValue(chr)
            val chrLength = chromosomeLengths.getValue(chr)
            val start = pileUpStart(record, chrLength, pileUpOptions.forwardShift, pileUpOptions.reverseShift)
            when (pileUpOptions.pileUpAlgorithm) {
                PileUpAlgorithm.START -> {
                    chrValues[start]++
                    sums[chr] = sums.getValue(chr) + 1
                }
                PileUpAlgorithm.MID_POINT -> {
                    val end = pileUpEnd(record, chrLength, pileUpOptions.forwardShift, pileUpOptions.reverseShift)
                    val midPoint = (start + end) / 2
                    chrValues[midPoint]++
                    sums[chr] = sums.getValue(chr) + 1
                }
                PileUpAlgorithm.LENGTH -> {
                    val end = pileUpEnd(record, chrLength, pileUpOptions.forwardShift, pileUpOptions.reverseShift)
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

    val pileUps = values.keys.map { chr ->
        chr to PileUp(
            chrLength = chromosomeLengths.getValue(chr),
            values = values.getValue(chr),
            sum = sums.getValue(chr)
        )
    }.toMap().toMutableMap()

    val pileUpSummary = pileUps.entries.joinToString("\n") {
        "chromosome ${it.key}: len=${it.value.chrLength} sum=${it.value.sum}"
    }
    log.info { "Pile-up complete with results: \n$pileUpSummary" }

    return pileUps
}

/**
 * Calculate the "start" value for our pile-up algorithm for the given record and settings.
 *
 * This start will actually be the record's end it's strand is minus
 */
private fun pileUpStart(record: SAMRecord, chrLength: Int, forwardShift: Int, reverseShift: Int): Int {
    return if (!record.readNegativeStrandFlag) {
        (record.start + forwardShift).withinBounds(0, chrLength)
    } else {
        (record.end + reverseShift).withinBounds(0, chrLength)
    }
}

/**
 * Calculate the "end" value for our pile-up algorithm for the given SAM record and settings.
 *
 * This end will actually be the record's start it's strand is minus
 */
private fun pileUpEnd(record: SAMRecord, chrLength: Int, forwardShift: Int, reverseShift: Int): Int {
    return if (!record.readNegativeStrandFlag) {
        (record.end + reverseShift).withinBounds(0, chrLength)
    } else {
        (record.start + forwardShift).withinBounds(0, chrLength)
    }
}

/**
 * @returns: if < start: start, if > end: end, else value
 */
private fun Int.withinBounds(start: Int, end: Int) =
    if (this < start) start else if (this > end) end else this